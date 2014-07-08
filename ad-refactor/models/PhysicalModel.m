classdef PhysicalModel
    properties
        % Unique identifier for the model
        name
        % The fluid model
        fluid
        % Operators used for construction of systems
        operators
        % Inf norm tolerance for nonlinear iterations
        nonlinearTolerance
        % Grid
        G
        
        % Maximum pressure change
        dpMax
        % Maximum saturation change
        dsMax
        
        % Water phase present
        water
        % Gas phase present
        gas
        % Oil phase present
        oil
        
        % Input data used to instansiate the model
        inputdata
    end
    
    methods
        function model = PhysicalModel(G, rock, fluid, varargin) %#ok
            model.dpMax = inf;
            model.dsMax = .2;
            model.nonlinearTolerance = 1e-6;
            model.inputdata = [];
            
            % Physical model
            model.G = G;
            model.fluid = fluid;
            
            model = merge_options(model, varargin{:});
            
            % Base class does not support any phases
            model.water = false;
            model.gas = false;
            model.oil = false;

        end
        
        function model = setupOperators(model, G, rock, varargin)
            % Set up divergence/gradient/transmissibility operators
            model.operators = setupSimComp(G, rock, varargin{:});
        end
        
        function [problem, state] = getEquations(model, state0, state, drivingForces, dt, varargin) %#ok
            % Get the equations governing the system
            error('Base class not meant for direct use')
        end
        
        function [state, report] = updateState(model, state, dx, drivingForces) %#ok
            % Update state based on non-linear increment
            error('Base class not meant for direct use')
        end
        
        function [convergence, values] = checkConvergence(model, problem, n)
            % Check convergence based on residual tolerances
            if nargin == 2
                n = inf;
            end
            
            values = norm(problem, n);
            convergence = all(values < model.nonlinearTolerance);
            
            if mrstVerbose()
                for i = 1:numel(values)
                    fprintf('%s (%s): %2.2e\t', problem.equationNames{i}, problem.types{i}, values(i));
                end
                fprintf('\n')
            end
        end
        
        function [state, report] = stepFunction(model, state, state0, dt, drivingForces, linsolve, nonlinsolve, varargin)
            % Make a single linearized timestep
            [problem, state] = model.getEquations(state0, state, dt, drivingForces, varargin{:});
            [convergence, values] = model.checkConvergence(problem);
            
            if ~convergence
                % Get increments for Newton solver
                [dx, ~, linearReport] = linsolve.solveLinearProblem(problem, model);
                
                % Let the non-linear solver decide what to do with the
                % increments to get the best convergence
                dx = nonlinsolve.stabilizeNewtonIncrements(problem, dx);
                
                % Finally update the state. The physical model knows which
                % properties are actually physical.
                [state, updateReport] = model.updateState(state, problem, dx, drivingForces);
                failure = any(cellfun(@(d) ~all(isfinite(d)), dx));
            else
                failure = false;
                [linearReport, updateReport] = deal(struct());
            end
            report = struct('LinearSolver', linearReport, ...
                            'UpdateState',  updateReport, ...
                            'Failure',      failure, ...
                            'Converged',    convergence, ...
                            'Residuals',    values);
        end
        
        function [gradient, result, report] = solveAdjoint(model, solver, getState,...
                                    getObjective, schedule, gradient, itNo)
            dt_steps = schedule.step.val;
            
            current = getState(itNo);
            before    = getState(itNo - 1);
            dt = dt_steps(itNo);
            
            lookupCtrl = @(step) schedule.control(schedule.step.control(step));
            [~, forces] = model.getDrivingForces(lookupCtrl(itNo));
            problem = model.getEquations(before, current, dt, forces, 'iteration', inf);
            
            if itNo < numel(dt_steps)
                after    = getState(itNo + 1);
                dt_next = dt_steps(itNo + 1);
                
                [~, forces_p] = model.getDrivingForces(lookupCtrl(itNo + 1));
                problem_p = model.getEquations(current, after, dt_next, forces_p,...
                                    'iteration', inf, 'reverseMode', true);
            else
                problem_p = [];
            end
            [gradient, result, rep] = solver.solveAdjointProblem(problem_p,...
                                        problem, gradient, getObjective(itNo), model);
            report = struct();
            report.Types = problem.types;
            report.LinearSolverReport = rep;
        end
        
        function [vararg, driving] = getDrivingForces(model, control) %#ok
            vararg = {};
            driving = struct('Wells', [], 'bc', [], 'src', []);
            
            if isfield(control, 'W') && ~isempty(control.W)
                vararg = [vararg, 'Wells', control.W];
                driving.Wells = control.W;
            end

            if isfield(control, 'bc') && ~isempty(control.bc)
                vararg = [vararg, 'bc', control.bc];
                driving.bc = control.bc;
            end
            
            if isfield(control, 'src') && ~isempty(control.src)
                vararg = [vararg, 'src', control.src];
                driving.src = control.src;
            end
        end
    end
end

