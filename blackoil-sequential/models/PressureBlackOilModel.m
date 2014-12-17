classdef PressureBlackOilModel < ThreePhaseBlackOilModel
    % Two phase oil/water system without dissolution
    properties

    end
    
    methods
        function model = PressureBlackOilModel(G, rock, fluid, varargin)
            
            model = model@ThreePhaseBlackOilModel(G, rock, fluid);

            model = merge_options(model, varargin{:});

            % Ensure simple tolerances
            model.useCNVConvergence = false;
        end
        
        function [problem, state] = getEquations(model, state0, state, dt, drivingForces, varargin)
            [problem, state] = pressureEquationBlackOil(state0, state, model,...
                            dt, ...
                            drivingForces,...
                            varargin{:});
            
        end
        function [convergence, values] = checkConvergence(model, problem, varargin)
            [convergence, values] = checkConvergence@PhysicalModel(model, problem, varargin{:});
            % Always make at least one update so that the problem actually changes.
            convergence = convergence && problem.iterationNo > 1;
        end
        
        function [state, report] = updateState(model, state, problem, dx, drivingForces)
            [state, report] = updateState@ReservoirModel(model, state, problem, dx, drivingForces);
            
            
            if 1
            state0 = state;
            state = updateWellCellSaturationsExplicit(model, state, problem, dx, drivingForces);
            
            
%             so = model.getProp(state, 'so');
%             sw = model.getProp(state, 'sw');
%             sg = model.getProp(state, 'sg');
% 
%             st = getCellStatusVO(state0, so, sw, sg, model.disgas, model.vapoil);
%             state = computeFlashBlackOil(state, state0, model, st);
            end
        end
        
        function [state, report] = updateAfterConvergence(model, state0, state, dt, drivingForces) %#ok
            state.s = state0.s;
            state.rv = state0.rv;
            state.rs = state0.rs;
            report = [];
        end
    end
end
