classdef WellModel < handle
%DEPRECATED Well model for three phase flow with black oil-style fluids
%
% SYNOPSIS:
%   wm = WellModel();
%
% DESCRIPTION:
%   This well model implements well equations, source terms and logic
%   related to controls switching for three phase black oil-like models
%   (possibly with dissolved gas and oil in the vapor phase).
%
% REQUIRED PARAMETERS:
%   None.
%
% OPTIONAL PARAMETERS:
%
%   (See documented class properties)
%
% RETURNS:
%   WellModel class instance.
%
% NOTE:
%   This class is deprecated!
%
% SEE ALSO:
%   ReservoirModel, ThreePhaseBlackOilModel

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

    properties
        % Index into list of phase pressures corresponding to the pressure
        % used for properties (typically oil pressure for black oil)
        referencePressureIndex
        % If enabled, injectors becoming producers and vice versa will not
        % be shut down.
        allowWellSignChange
        % Allow different sign per completion rate
        allowCrossflow
        % Allow wells to switch controls if limits are hit (otherwise shut
        % down).
        allowControlSwitching
        % Extra output to command line
        verbose
        % Include extra fields in well sol (gas / oil ratio, rates at
        % reservoir conditions etc).
        detailedOutput
        
        %  Properties used during calculations. These are not set by the
        %  constructor, but rather used during function calls. All values
        %  are given per well cell.
        
        % Cell array of number of phases b-factors
        bfactors
        % Surface density for each phase
        surfaceDensities
        % Cell mobilities for each phase
        mobilities
        % Phase pressures
        pressure
        % The pressure used as the global pressure (typically oil pressure)
        referencePressure
        % max Rs / maxRv, if applicable
        maxComponents
        % Rs and Rv in cells
        components
        % Saturations per phase in well cells
        saturations
        % Indicator for nonlinear iteration number. Some features will only
        % be computed during the first nonlinear iteration of a timestep.
        nonlinearIteration
        % The current well controls
        W
        % Mapping from perforation index to well index. perf2well(ic) returns the well index of the
        % perforation ix
        perf2well
        % Inverse mapping for perf2wll. Rw(:, iw) returns the perforations that belongs to the well iw, in
        % form of a logical vector
        Rw
    end
    
    methods
        function wmodel = WellModel(varargin)
            wmodel.allowWellSignChange   = false;
            wmodel.allowCrossflow        = true;
            wmodel.allowControlSwitching = true;
            wmodel = merge_options(wmodel, varargin{:});
        end
        
        function [sources, wellEqs, controlEqs, wc, wellSol, sources_reservoir] = ...
                            computeWellFlux(wellmodel, model,...
                                            W, wellSol, bhp,...
                                            currentFluxes, pressure,...
                                            surfaceDensities, bfactors, ...
                                            mob, satvals, ...
                                            compvals, varargin)
                                        
            
            wellmodel.verbose = mrstVerbose();
            wellmodel.maxComponents = {};
            wellmodel.nonlinearIteration = nan;
            wellmodel.referencePressureIndex = 2;
            wellmodel.detailedOutput        = model.extraWellSolOutput;
            
            wellmodel = merge_options(wellmodel, varargin{:});
            
            % Store all current variables in the current well model
            % instance temporarily
            wellmodel.bfactors = bfactors;
            wellmodel.surfaceDensities = surfaceDensities;
            wellmodel.mobilities = mob;
            wellmodel.saturations = satvals;
            wellmodel.components = compvals;
            wellmodel.W = W;
            [wellmodel.perf2well, wellmodel.Rw] = getPerfToWellMap(wellmodel);
            
            if isempty(W)
                sources = {};
                controlEqs = {};
                return
            end
            nsat = numel(model.saturationVarNames);
            
            if ~iscell(pressure)
                % Support single reference pressure
                p = cell(nsat, 1);
                [p{:}]  = deal(pressure);
            else
                p = pressure;
            end
            wellmodel.referencePressure = p{wellmodel.referencePressureIndex};
            wellmodel.pressure = p;
            clear pressure
            
            % Assertions to indirectly document the implementation
            assert(numel(p) == nsat, ...
                'Number of saturation pressure did not match number of saturation variables!');
            assert(numel(mob) == nsat,...
                'Number of saturation mobilities did not match number of saturation variables!');
            assert(numel(currentFluxes) == nsat, ...
                'Number of phase fluxes did not match number of saturation variables!');
            assert(numel(satvals) == nsat,...
                'Number of provided saturation values did not match number of saturation variables!');
            assert(numel(surfaceDensities) == nsat && numel(bfactors) == nsat, ...
                'Number of densities or bfactors did not match number of present phases!');
            
            % Update well pressure
            [wellSol, currentFluxes, bhp] = wellmodel.updatePressure(wellSol, currentFluxes, bhp, model);
            % Update well limits
            [wellSol, currentFluxes, bhp] = wellmodel.updateLimits(wellSol, currentFluxes, bhp, model);
            
            if false
                wellmodel = setUpstreamMobility(wellmodel, model, wellSol, bhp, currentFluxes);
            end
            
            % Set up the actual equations
            [wellEqs, controlEqs, sources, sources_reservoir, wellSol] =...
                wellmodel.assembleEquations(wellSol, currentFluxes, bhp, model);
            
            if wellmodel.detailedOutput
                wellSol = wellmodel.updateWellSolStatistics(wellSol, sources, model);
            end
            % Check for multiple perforations in same cells to account for
            % a shortcoming in MATLABs indexing behavior (repeated indices
            % are not summed, they are overwritten).
            wc = vertcat(W.cells);
            [wc, sources] = wellmodel.handleRepeatedPerforatedcells(wc, sources);
        end
        
        function [wellSol, q_s, bhp] = updateLimits(wellmodel, wellSol, q_s, bhp, model)
            % Update the switched well controls based on current values
            if ~wellmodel.allowControlSwitching
                return
            end
            % First pad all values so that the well and fluxes are three
            % phase in appearance (with zeros for non-present values)
            [q_s, W_3ph, active] = padRatesAndCompi(q_s, wellmodel.W, model);
            
            [wellSol, withinLims] = updateSwitchedWellControls(wellmodel, model, wellSol, bhp, q_s);
            
            if all(withinLims)
                q_s = q_s(active);
                return
            end
            
            for k = 1:numel(wellSol)
                w = W_3ph(k);
                if ~withinLims(k)
                    v  = wellSol(k).val;
                    switch wellSol(k).type
                        case 'bhp'
                            bhp = assignValue(bhp, v, k);
                        case 'rate'
                            q_s{1} = assignValue(q_s{1}, v*w.compi(1), k);
                            q_s{2} = assignValue(q_s{2}, v*w.compi(2), k);
                            if numel(q_s)>2
                                q_s{3} = assignValue(q_s{3}, v*w.compi(3), k);
                            end
                        case 'orat'
                            q_s{2} = assignValue(q_s{2}, v, k);
                        case 'wrat'
                            q_s{1} = assignValue(q_s{1}, v, k);
                        case 'grat'
                            q_s{3} = assignValue(q_s{3}, v, k);
                    end % No good guess for qOs, etc...
                end
            end
            
            q_s = q_s(active);
        end
        
        function [wellSol, flux, bhp] = updatePressure(wellmodel, wellSol, flux, bhp, model)
            % Update pressure drop along well bore. Will typically only do
            % work at the first nonlinear iteration.
            if isnan(wellmodel.nonlinearIteration) || wellmodel.nonlinearIteration < 0
                warning(['Iteration number is not passed on to well model,', ...
                         'this may indicate welbore pressure-drop will never be updated']);
            end
            
            if wellmodel.nonlinearIteration == 1
                wellSol = updateConnectionDP(wellmodel, model, wellSol);
            end
            
            if all(isfinite(double(bhp)))
                return
            end
            % Otherwise we have some initialization to do
            pd   = double(wellmodel.referencePressure);
            for k = 1:numel(wellSol)
                if ~isfinite(wellSol(k).bhp)
                    v = pd(k) - wellSol(k).cdp(1);
                    if isempty(wellSol(k).sign)
                        dispif(wellmodel.verbose, 'No sign found in wellsol...\n');
                    else
                        v = v + 5*wellSol(k).sign*barsa;
                    end
                    bhp = assignValue(bhp, v, k);
                    wellSol(k).bhp = v;
                end
            end
        end
        
        
        function  [wellEqs, controlEqs, cq_s, cq_r, wellSol] = assembleEquations(wellmodel,...
                                                wellSol, q_s, bhp, model)
            % Assemble well model equations (Peaceman type), well control
            % equations, completion surface and reservoir rates, plus an
            % updated wellSol corresponding to the current limits.
            [wellEqs, cq_s, mix_s, status, cstatus, cq_r] = ...
                            computeWellContributionsNew(wellmodel, model, wellSol, bhp, q_s);
            controlEqs =  setupWellControlEquations(wellSol, bhp, q_s, status, mix_s, model);
            
            
            % Update well properties which are not primary variables
            toDouble = @(x)cellfun(@double, x, 'UniformOutput', false);
            cq_sDb = cell2mat(toDouble(cq_s));
            
            for wnr = 1:numel(wellSol)
                ix = wellmodel.perf2well == wnr;
                wellSol(wnr).cqs     = cq_sDb(ix,:);
                wellSol(wnr).cstatus = cstatus(ix);
                wellSol(wnr).status  = status(wnr);
            end
        end
        
        function ws = updateWellSolStatistics(wellmodel, ws, sources, model)
            % Store extra output, typically black oil-like
            p2w = wellmodel.getPerfToWellMap();
            
            gind = model.getPhaseIndex('G');
            oind = model.getPhaseIndex('O');
            wind = model.getPhaseIndex('W');
            if ~isempty(wellmodel.bfactors)
                bf  = cellfun(@double, wellmodel.bfactors, 'UniformOutput', false);
            else
                bf = cell(numel(sources), 1);
                [bf{:}] = deal(ones(size(p2w)));
            end
            src = cellfun(@double, sources, 'UniformOutput', false);
            for i = 1:numel(ws)
                % Store reservoir fluxes and total fluxes
                ws(i).qTs = 0;
                ws(i).qTr = 0;
                if model.gas
                    tmp = sum(src{gind}(p2w == i)./bf{gind}(p2w == i));
                    ws(i).qGr = tmp;
                    ws(i).qTr = ws(i).qTr + tmp;
                    ws(i).qTs = ws(i).qTs + ws(i).qGs;
                end
                
                if model.oil
                    tmp = sum(src{oind}(p2w == i)./bf{oind}(p2w == i));
                    ws(i).qOr = tmp;
                    ws(i).qTr = ws(i).qTr + tmp;
                    ws(i).qTs = ws(i).qTs + ws(i).qOs;
                end
                
                if model.water
                    tmp = sum(src{wind}(p2w == i)./bf{wind}(p2w == i));
                    ws(i).qWr = tmp;
                    ws(i).qTr = ws(i).qTr + tmp;
                    ws(i).qTs = ws(i).qTs + ws(i).qWs;
                end
                
                % Phase cuts - fraction of reservoir conditions
                if model.water && model.oil
                    ws(i).wcut = ws(i).qWs./(ws(i).qWs + ws(i).qOs);
                end
                
                if model.gas
                    ws(i).gcut = ws(i).qGs./ws(i).qTs;
                end
                
                if model.oil
                    ws(i).ocut = ws(i).qOs./ws(i).qTs;
                end
                
                % Gas/oil ratio
                if model.gas && model.oil
                    ws(i).gor = ws(i).qGs/ws(i).qOs;
                end
            end
        end
        
        function wellmodel = setUpstreamMobility(wellmodel, model, wellSol, bhp, q_s)
            % Upstream weighting of injection mobilities
            wellStatus = vertcat(wellmodel.W.status);
            % Well total volume rate at std conds:
            qt_s = q_s{1};
            nph = numel(q_s);
            for ph = 2:nph
                qt_s = (qt_s + q_s{ph}).*wellStatus;
            end
            inj = double(qt_s) > 0;
            p2w = wellmodel.getPerfToWellMap();            
            compi = vertcat(wellmodel.W.compi);
            perfcompi = compi(p2w, :);

            RMw    = sparse((1:numel(p2w))', p2w, 1, numel(p2w), numel(wellmodel.W));
            drawdown = -(RMw*bhp+vertcat(wellSol.cdp)) + wellmodel.referencePressure;
            
            
            [kr, mu, sat] = deal(cell(1, nph));
            for i = 1:nph
                sat{i} = perfcompi(:, i);
            end
            [kr{:}] = model.evaluateRelPerm(sat);
            
            
            f = model.fluid;
            ix = 1;
            if model.water
                mu{ix} = f.muW(drawdown);
                ix = ix + 1;
            end
            
            if model.oil
                if isprop(model, 'disgas') && model.disgas
                    isgas = sat{ix} > 0;
                    rs = model.fluid.rsSat(drawdown);
                    mu{ix} = f.muO(drawdown, rs.*isgas, isgas);
                else
                    if isfield(f, 'BOxmuO')
                        mu{ix} = f.BOxmuO(drawdown).*f.bO(drawdown);
                    else
                        mu{ix} = f.muO(drawdown);
                    end
                end
                ix = ix + 1;
            end
            
            if model.gas
                if isprop(model, 'vapoil') && model.vapoil
                    % Strange case
                    isoil = sat{ix} > 0;
                    rv = model.fluid.rvSat(drawdown);
                    mu{ix} = f.muG(drawdown, rv.*isoil, isoil);
                else
                    mu{ix} = f.muG(drawdown);
                end
            end
            mob = cellfun(@(x, y) x./y, kr, mu, 'UniformOutput', false);
            
            for i = 1:nph
                injperf = inj(p2w);
                wellmodel.mobilities{i}(injperf) = mob{i}(injperf).*perfcompi(injperf, i);
            end
        end
        
        function varargout = getPerfToWellMap(wellmodel)
            % Outsource this, but it could be overriden
            varargout = cell(nargout, 1);
            [varargout{:}] = getPerforationToWellMapping(wellmodel.W);
        end
    end
    methods (Static)
        function [wc, cqs] = handleRepeatedPerforatedcells(wc, cqs)
            % This function treats repeated indices in wc (typically due to
            % multiple wells intersecting a single cell). The output will
            % have no repeats in wc, and add together any terms in cqs.
            [c, ic, ic] = uniqueStable(wc);                     %#ok<ASGLU>
            if numel(c) ~= numel(wc)
                A = sparse(ic, (1:numel(wc))', 1, numel(c), numel(wc));
                wc = c;
                for k=1:numel(cqs)
                    cqs{k} = A*cqs{k};
                end
            end
        end
        
        function [rw, rsatw] = getResSatWell(model, cells, rs, rv, rsSat, rvSat)
            nperf = numel(cells);
            if model.disgas
                % Sample the cells
                rsw = rs(cells); 
                rsSatw = rsSat(cells);
            else
               % rs supposed to be scalar in this case
                rsw = ones(nperf,1)*mean(rs); 
                rsSatw = ones(nperf,1)*mean(rsSat); 
            end
            if model.vapoil
                rvw = rv(cells); 
                rvSatw = rvSat(cells);
            else
                % rv supposed to be scalar in this case
                rvw = ones(nperf,1)*mean(rv); 
                rvSatw = ones(nperf,1)*mean(rvSat); 
            end
            
            rw = {rsw, rvw};
            rsatw = {rsSatw, rvSatw};
        end
        
        function [eqs, names, types] = createReverseModeWellEquations(model, wellSol, sampleVariable)
            nph = sum(model.getActivePhases());
            nw = numel(wellSol);
            zw = double2ADI(zeros(nw,1), sampleVariable);
            
            [eqs, names, types] = deal(cell(1, nph+1));
            [eqs{:}]   = deal(zw);
            [names{:}] = deal('empty');
            [types{:}] = deal('none');
        end
    end
end



