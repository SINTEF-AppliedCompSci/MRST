classdef TransportOverallCompositionModel < OverallCompositionCompositionalModel
    % Two phase oil/water system without dissolution
    properties
        staticUpwind
        upwindType
    end
    
    methods
        function model = TransportOverallCompositionModel(G, rock, fluid, compfluid, varargin)
            
            model = model@OverallCompositionCompositionalModel(G, rock, fluid, compfluid);
            
            model.staticUpwind = false;
            model.upwindType  = 'potential';
            model.useIncTolComposition = false;

            model = merge_options(model, varargin{:});
        end
        
        function [problem, state] = getEquations(model, state0, state, dt, drivingForces, varargin)
            if ~isfield(state, 'sT')
                state.sT = sum(state.s, 2);
            end
            if ~isfield(state0, 'sT')
                state0.sT = sum(state0.s, 2);
            end
            [problem, state] = transportEquationOverallComposition(state0, state, model,...
                            dt, ...
                            drivingForces,...
                            varargin{:});
            
        end
        
        function [state, report] = updateState(model, state, problem, dx, drivingForces)
            isS = strcmpi(problem.primaryVariables, 'sT');
            ds = dx{isS};
            state = model.updateStateFromIncrement(state, ds, problem, 'sT', inf, inf);
            problem.primaryVariables = problem.primaryVariables(~isS);
            dx = dx(~isS);
            
            [state, report] = updateState@OverallCompositionCompositionalModel(model, state, problem, dx, drivingForces);
        end
        function state = validateState(model, state)
            state.sT = ones(model.G.cells.num, 1);
        end
        
        function [fn, index] = getVariableField(model, name, varargin)
            switch(lower(name))
                case {'st'}
                    fn = 'sT';
                    index = 1;
                otherwise
                    [fn, index] = getVariableField@OverallCompositionCompositionalModel(model, name, varargin{:});
            end
        end
        
        function [state, report] = updateAfterConvergence(model, state0, state, dt, drivingForces)
            [state, report] = updateAfterConvergence@OverallCompositionCompositionalModel(model, state0, state, dt, drivingForces);
            state.s = bsxfun(@times, state.s, state.sT);
            state = rmfield(state, 'sT');
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
