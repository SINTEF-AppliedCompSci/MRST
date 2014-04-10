classdef physicalModel
    properties
        name
        fluid
        operators
        nonlinearTolerance
        G
        
        % capping defaults
        dpMax
        dsMax
    end
    
    methods
        function model = physicalModel()
            model.dpMax = inf;
            model.dsMax = .2;
            model.nonlinearTolerance = 1e-8;
        end
        
        function model = setupOperators(model, G, rock, varargin)
            model.operators = setupSimComp(G, rock, varargin{:});
        end
        
        function eqs = getEquations(model, state0, state, drivingForces, dt, varargin) %#ok
            error('Base class not meant for direct use')
        end
        
        function state = updateState(model, state, dx, drivingForces) %#ok
            error('Base class not meant for direct use')
        end
        
        function [convergence, values] = checkConvergence(model, problem, n)
            % Basic convergence check
            if nargin == 2
                n = inf;
            end
            
            values = norm(problem, n);
            convergence = all(values < model.nonlinearTolerance);
            
            if mrstVerbose()
                for i = 1:numel(values)
                    fprintf('%s (%s): %2.2e\t', problem.names{i}, problem.types{i}, values(i));
                end
                fprintf('\n')
            end
        end
        
        function [state, convergence] = stepFunction(model, state, state0, dt, drivingForces, solver, varargin)
            problem = model.getEquations(state0, state, dt, drivingForces, varargin{:});
            convergence = model.checkConvergence(problem);
            if ~convergence
                dx = solver.solveLinearProblem(problem);
                state = model.updateState(state, dx, drivingForces); %%
            end
        end
    end
end