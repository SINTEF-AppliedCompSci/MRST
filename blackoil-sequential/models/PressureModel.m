classdef PressureModel < WrapperModel
    properties
        incTolPressure = 1e-3;
        useIncTol = true;
        reductionStrategy = 'analytic';
        pressureIncTolType = 'relative';
    end
    
    methods
        function model = PressureModel(parent, varargin)
            if isprop(parent, 'useCNVConvergence')
                parent.useCNVConvergence = false;
            end
            model = model@WrapperModel(parent);
            model = merge_options(model, varargin{:});
            model.AutoDiffBackend = parent.AutoDiffBackend;
        end
        
        function [state, names, origin] = getStateAD(model, state, init)
            if nargin < 3
                init = true;
            end
            % Get the AD state for this model
            [vars, names, origin] = model.getPrimaryVariables(state);
            switch model.reductionStrategy
                case 'analytic'
                    isP = strcmp(names, 'pressure');
                    origP = origin{isP};
                    keep = isP | cellfun(@(x) ~strcmp(x, origP), origin);
                case 'numerical'
                    
            end
            if init
                [vars{keep}] = model.AutoDiffBackend.initVariablesAD(vars{keep});
            end
            state = model.initStateAD(state, vars, names, origin);
            % Not all were AD-initialized
            names = names(keep);
            origin = origin(keep);
        end

        function [eqs, names, types, state] = getModelEquations(pmodel, state0, state, dt, drivingForces)
            model = pmodel.parentModel;
            
            [ceqs, flux] = model.FluxDiscretization.componentConservationEquations(model, state, state0, dt);
            w = model.getProp(state, 'PressureReductionFactors');

            src = model.FacilityModel.getComponentSources(state);
            % Assemble equations and add in sources
            pressure_equation = 0;
            for i = 1:numel(ceqs)
                if ~isempty(src.cells)
                    ceqs{i}(src.cells) = ceqs{i}(src.cells) - src.value{i};
                end
                pressure_equation = pressure_equation + w{i}.*model.operators.AccDiv(ceqs{i}, flux{i});
            end
            % Get facility equations
            [weqs, wnames, wtypes, state] = model.FacilityModel.getModelEquations(state0, state, dt, drivingForces);
            eqs = [{pressure_equation}, weqs];
            names = ['pressure', wnames];
            types = ['cell', wtypes];
        end

        
        function [problem, state] = getEquations(model, state0, state, dt, drivingForces, varargin)
            [problem, state] = getEquations@PhysicalModel(model, state0, state, dt, drivingForces, varargin{:});
        end
        
        function [state, report] = updateState(model, state, problem, dx, drivingForces)
            p0 = state.pressure;
            [state, report] = model.parentModel.updateState(state, problem, dx, drivingForces);
            switch lower(model.pressureIncTolType)
                case 'relative'
                    range = max(p0) - min(p0);
                    if range == 0
                        range = 1;
                    end
                    state.dpRel = (state.pressure - p0)./range;
                case 'absolute'
                    state.dpRel = state.pressure - p0;
                case 'norm'
                    state.dpRel = (state.pressure - p0)/norm(state.pressure, inf);
            end
        end
        
        function  [convergence, values, names] = checkConvergence(model, problem)
            [convergence, values, names] = model.parentModel.checkConvergence(problem);
            if ~isnan(problem.iterationNo) && model.useIncTol
                if problem.iterationNo  > 1
                    values(1) = norm(problem.state.dpRel, inf);
                else
                    values(1) = inf;
                end
                convergence = [values(1) < model.incTolPressure, convergence(2:end)];
                names{1} = 'Delta P';
            end
        end

        function [state, report] = updateAfterConvergence(model, varargin)
            [state, report] = updateAfterConvergence@WrapperModel(model, varargin{:});
            if isfield(state, 'statePressure')
                state = rmfield(state, 'statePressure');
            end
            if isfield(state, 'sT')
                state.sT = sum(value(state.s), 2);
            end
            state.statePressure = state;
        end
        
        function rhoS = getSurfaceDensities(model)
            rhoS = model.parentModel.getSurfaceDensities();
        end
    end
end

%{
Copyright 2009-2019 SINTEF Digital, Mathematics & Cybernetics.

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
