classdef FacilityModelSolvent < FacilityModel
    
    methods
        function model = FacilityModelSolvent(reservoirModel, varargin)
            model = model@FacilityModel(reservoirModel, varargin{:});
        end
        
        function names = getBasicPrimaryVariableNames(model)
            % Basic primary variables are phase rates + bhp for active
            % phases in the model.
            actPh = model.ReservoirModel.getActivePhases();
            names = {'qWs', 'qOs', 'qGs', 'qSs', 'bhp'};
            names = names([actPh, true]);
        end

        function [rates, bhp, names] = getBasicPrimaryVariables(model, wellSol)
            % Get phase rates + bhp for active phases
            if model.getNumberOfWells() == 0
                [rates, names] = deal({});
                bhp = [];
            else
                actPh = model.ReservoirModel.getActivePhases();
                bhp = vertcat(wellSol.bhp);
                qWs = vertcat(wellSol.qWs);
                qOs = vertcat(wellSol.qOs);
                qGs = vertcat(wellSol.qGs);
                qSs = vertcat(wellSol.qSs);
                rates = {qWs, qOs, qGs, qSs};
                rates = rates(actPh);

                names = model.getBasicPrimaryVariableNames();
            end
        end
        
        function model = setupWells(model, W, wellmodels)
            % Set up wells for changed controls or first simulation step
            nw = numel(W);
            if model.getNumberOfWells == 0
                % First time setup
                [pvars, eqnames, eqtypes] = deal(cell(nw, 1));
                model.WellModels = cell(nw, 1);
                for i = 1:nw
                    % Set up models. SimpleWell for the time being
                    if nargin < 3
                        if isfield(W(i), 'isMS') && W(i).isMS
                            wm = MultisegmentWell(W(i));
                        else
                            wm = SimpleWellSolvent(W(i));
                        end
                    else
                        wm = wellmodels{i};
                    end
                    wm.dsMaxAbs = model.ReservoirModel.dsMaxAbs;
                    wm.dpMaxAbs = model.ReservoirModel.dpMaxAbs;
                    wm.dpMaxRel = model.ReservoirModel.dpMaxRel;
                    
                    if isfield(W(i), 'vfp_index')
                        vfp_ix = W(i).vfp_index;
                        if vfp_ix > 0
                            if wm.isInjector()
                                vfp = model.VFPTablesInjector{vfp_ix};
                            else
                                vfp = model.VFPTablesProducer{vfp_ix};
                            end
                            wm.VFPTable = vfp;
                        end
                    end
                    % Get the added primary variables for this well, plus
                    % the equations and equation types it adds
                    pvars{i} = wm.getExtraPrimaryVariableNames(model.ReservoirModel);
                    [eqnames{i}, eqtypes{i}] = wm.getExtraEquationNames(model.ReservoirModel);
                    model.WellModels{i} = wm;
                end
                % Combine the different equations and types added by the
                % different wells into a canonical ordering.
                model.addedPrimaryVarNames = uniqueStable([pvars{:}]);
                [model.addedEquationNames, keepix] = uniqueStable([eqnames{:}]);

                etypes = [eqtypes{:}];
                model.addedEquationTypes = etypes(keepix);
            else
                assert(model.getNumberOfWells == nw, ...
                    'Number of wells in facility model has changed during simulation')
                for i = 1:nw
                    % Update with new wells. Typically just a quick
                    % overwrite of existing wells
                    model.WellModels{i} = model.WellModels{i}.updateWell(W(i));
                end
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
