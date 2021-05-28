classdef UniformFacilityModel < FacilityModel
    % Simplified facility model which is sometimes faster
    properties
        
    end
    
    methods
        
        function model = UniformFacilityModel(varargin)
            model = model@FacilityModel(varargin{:});
        end
        
        function [sources, wellSystem, wellSol] = getWellContributions(model, wellSol0, wellSol, wellvars, wellMap, p, mob, rho, dissolved, comp, dt, iteration)
            if isnan(iteration) || iteration < 0
                warning(['Iteration number is not passed on to well model,', ...
                         'this may indicate wellbore pressure-drop will never be updated']);
            end
            actWellIx = model.getIndicesOfActiveWells(wellSol);
            nw = numel(actWellIx);

            % Get the additional equations not implemented in the minimal
            % "SimpleWell" class.
            enames = model.addedEquationNames;
            cnames = model.ReservoirModel.getComponentNames();
            ncomp = numel(cnames);
            assert(ncomp == 0, 'UniformFacilityModel does not support components.');
            
            n_extra = numel(enames);
            assert(n_extra == 0);
            resModel = model.ReservoirModel;

            addedVars = model.addedPrimaryVarNames;
            varmaps = cell(1, numel(addedVars));
            for varNo = 1:numel(addedVars)
                varmaps{varNo} = model.getWellVariableMap(addedVars{varNo}, wellSol);
            end
            
            isBH = wellMap.isBHP;
            isQ = wellMap.isRate;
            emap = wellMap.extraMap;
            
            bhp = wellvars{isBH};
            qWell = wellvars(isQ);
            wellvars = wellvars(~(isBH | isQ));
            
            wc = model.getWellCells(actWellIx);
            [basenames, basetypes] = model.WellModels{1}.getWellEquationNames(resModel);
            
            toDouble = @(x)cellfun(@value, x, 'UniformOutput', false);
            
            p_d = value(p);
            mob_d = toDouble(mob);
            rho_d = toDouble(rho);
            wellvars_d = toDouble(wellvars);
            dissolved_d = cellfun(toDouble, dissolved, 'UniformOutput', false);
            for i = 1:nw
                wellNo = actWellIx(i);
                wm = model.WellModels{wellNo};
                ws = wellSol(wellNo);
                ws0 = wellSol0(wellNo);

                W = wm.W;
                packed = packPerforationProperties(W, p_d, mob_d, rho_d, dissolved_d, comp, wellvars_d, addedVars, varmaps, emap, i);
                qw = cellfun(@(x) x(i), qWell, 'uniformoutput', false);
                bh = bhp(i);
                % Update pressure
                ws = wm.updateConnectionPressureDrop(ws0, ws, resModel, qw, bh, packed, dt, iteration);
                % Update limits
                [qw, bh, ws, ok] = wm.updateLimits(ws0, ws, resModel, qw, bh, packed, dt, iteration);
                if ~ok
                    bhp(i) = bh;
                    for phNo = 1:numel(qw)
                        qWell{phNo}(i) = qw{phNo};
                    end
                end
                wellSol(wellNo) = ws;
            end
            rhoS = model.ReservoirModel.getSurfaceDensities();
            extract = @(V) cellfun(@(x) x(wc), V, 'UniformOutput', false);

            p_w = p(wc);
            rho_w = extract(rho);
            mob_w = extract(mob);
            dissolved_w = dissolved;
            for i = 1:numel(dissolved_w)
                for j = 1:numel(dissolved_w{i})
                    if ~isempty(dissolved_w{i}{j})
                        dissolved_w{i}{j} = dissolved_w{i}{j}(wc);
                    end
                end
            end
            
            b_w = phaseDensitiesTobfactor(rho_w, rhoS, dissolved_w);
            if isa(model.ReservoirModel, 'ThreePhaseBlackOilModel') && ...
                  (model.ReservoirModel.disgas || model.ReservoirModel.vapoil)
                w = model.ReservoirModel.water;
                % RS, then RV
                rs = dissolved_w{2 + w}{1 + w};
                rv = dissolved_w{1 + w}{2 + w};
                
                if isempty(rs)
                    rs = 0;
                end
                if isempty(rv)
                    rv = 0;
                end
                dissolved_w = {rs, rv};
            end
            
            W = model.getWellStruct(actWellIx);
            
            wellmodel = struct('W', W ,'referencePressure', p_w,...
                        'bfactors', {b_w}, 'components', {dissolved_w}, ...
                        'mobilities', {mob_w});
            [wellmodel.perf2well, wellmodel.Rw] = getPerforationToWellMapping(W);
            wellmodel.allowWellSignChange = model.WellModels{1}.allowSignChange;
            wellmodel.allowCrossflow = model.WellModels{1}.allowCrossflow;
            wellmodel.water = model.ReservoirModel.water;
            
            
            [eqs, srcSurfCond, mix_s, status, cstatus, srcVol] = ...
                            computeWellContributionsNew(wellmodel, model.ReservoirModel, wellSol(actWellIx), bhp, qWell);
            ctrleq =  setupWellControlEquations(wellSol(actWellIx), bhp, qWell, status, mix_s, model.ReservoirModel);
            srcMass = srcSurfCond;
            for i = 1:numel(srcMass)
                srcMass{i} = srcMass{i}*rhoS(i);
            end
            
            names = basenames;
            types = basetypes;
            [wc, srcMass, srcVol] = model.handleRepeatedPerforatedcells(wc, srcMass, srcVol);
            wellSystem = struct('wellEquations', {eqs}, ...
                                'names',  {names}, ...
                                'types', {types}, ...
                                'controlEquation', ctrleq);
            sources = struct('phaseMass',   {srcMass}, ...
                             'phaseVolume', {srcVol}, ...
                             'components',  {{}}, ...
                             'sourceCells', wc);
            if model.ReservoirModel.extraWellSolOutput
                wellSol = model.setWellSolStatistics(wellSol, sources);
            end
            
            cq_sDb = cell2mat(toDouble(srcSurfCond));

            for i = 1:numel(actWellIx)
                wnr = actWellIx(i);
                ix = wellmodel.perf2well == i;
                wellSol(wnr).cqs     = cq_sDb(ix,:);
                wellSol(wnr).cstatus = cstatus(ix);
                wellSol(wnr).status  = status(i);
            end
        end
    end
end

%{
Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}
