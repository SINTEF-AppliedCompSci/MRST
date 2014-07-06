classdef nonlinearSolver < handle
    properties
        maxIterations
        maxSubsteps
        linearSolver
        verbose
        isAdjoint
        timeStepSelector
        
        % Related to stabilization of newton steps
        relaxationParameter
        relaxationType
        previousIncrement
        useRelaxation
    end
    
    methods
        function solver = nonlinearSolver(varargin)            
            solver.maxIterations = 25;
            solver.verbose       = mrstVerbose();
            solver.maxSubsteps   = 32;
            solver.isAdjoint     = false;
            solver.linearSolver  = [];
            
            solver.relaxationParameter = 1;
            solver.relaxationType = 'dampen';
            solver.useRelaxation = false;
            
            solver = merge_options(solver, varargin{:});
            
            if isempty(solver.linearSolver)
                solver.linearSolver = mldivideSolverAD();
            end
            
            if isempty(solver.timeStepSelector)
                solver.timeStepSelector = SimpleTimeStepSelector();
            end
        end
        
        function [state, report] = solveTimestep(solver, state0, dT, model, varargin)
            % Solve a timestep for a non-linear system using one or more substeps
            drivingForces = struct('Wells', [],...
                                   'bc',    [],...
                                   'src',   []);
            
            drivingForces = merge_options(drivingForces, varargin{:});
            
            
            converged = false;
            done = false;
            
            dt = dT;
            
            % Number of nonlinear iterations total
            itCount = 0;
            % Number of ministeps due to cutting 
            ministepNo = 1;
            % Number of steps
            stepCount = 0;
            % Number of accepted steps
            acceptCount = 0;
            
            t_local = 0;
            
            isFinalMinistep = false;
            state0_inner = state0;
            
            reports = {};
            
            state = state0;
            
            % Let the step selector know that we are at start of timestep
            % and what the current driving forces are
            stepsel = solver.timeStepSelector;
            stepsel.newControlStep(drivingForces);
            while ~done
                dt = stepsel.pickTimestep(dt, model, solver);
                
                if t_local + dt >= dT
                    % Ensure that we hit report time
                    isFinalMinistep = true;
                    dt = dT - t_local;
                end
                [state, nonlinearReports, converged, its] = ...
                        solveMinistep(solver, model, state, state0_inner, dt, drivingForces);
                
                % Store timestep info
                clear tmp;
                tmp.NonlinearReport = nonlinearReports;
                tmp.LocalTime = t_local + dt; 
                tmp.Converged = converged;
                tmp.Timestep = dt;
                tmp.Iterations = its;
                
                reports{end+1} = tmp; %#ok 
                stepsel.storeTimestep(tmp);
                
                % Keep total itcount so we know how much time we are
                % wasting
                itCount = itCount + its;
                stepCount = stepCount + 1;
                if converged
                    t_local = t_local + dt;
                    state0_inner = state;
                    acceptCount = acceptCount + 1;
                else
                    state = state0_inner;
                    % Beat timestep with a hammer
                    warning('Solver did not converge, cutting timestep')
                    ministepNo = 2*ministepNo;
                    dt = dt/2;
                    if ministepNo > solver.maxSubsteps
                        error('Did not find a solution. Reached maximum amount of substeps');
                    end
                    isFinalMinistep = false;
                end
                done = isFinalMinistep && converged;
            end
            dispif(solver.verbose, ...
                'Solved timestep with %d accepted ministeps (%d rejected, %d total iterations)\n',...
                acceptCount, stepCount - acceptCount, itCount);
            
            % Truncate reports from step functions
            reports = reports(~cellfun(@isempty, reports));
            report = struct('Iterations',       itCount,...
                            'Converged',        converged,...
                            'MinistepCount', 	ministepNo);
            % Add seperately because struct constructor interprets cell
            % arrays as repeated structs.
            report.StepReports = reports;
        end
        
        function dx = stabilizeNewtonIncrements(solver, problem, dx)
            dx_prev = solver.previousIncrement;
            
            w = solver.relaxationParameter;
            if w == 1
                return
            end
            
            switch(lower(solver.relaxationType))
                case 'dampen'
                    for i = 1:numel(dx)
                        dx{i} = dx{i}*w;
                    end
                case 'sor'
                    if isempty(dx_prev)
                        return
                    end
                    for i = 1:numel(dx)
                        dx{i} = dx{i}*w + (1-w)*dx_prev{i};
                    end
                case 'none'
                    
                otherwise
                    error('Unknown relaxationType: Valid options are ''dampen'', ''none'' or ''sor''');
            end
            solver.previousIncrement = dx;
        end
    end
end

function [state, reports, converged, its] = solveMinistep(solver, model, state, state0, dt, drivingForces)
    reports = cell(solver.maxIterations, 1);
    omega0 = solver.relaxationParameter;
    
    r = [];
    for i = 1:solver.maxIterations
        [state, stepReport] = ...
            model.stepFunction(state, state0, dt, drivingForces, ...
                               solver.linearSolver, solver, 'iteration', i);
        converged  = stepReport.Converged;
        if converged
            break
        end
        reports{i} = stepReport;
        
        
        if i > 1 && solver.useRelaxation
            % Check relative improvement
            improvement = (r_prev - stepReport.Residuals)./r_prev;
            % Remove NaN / inf due to converged / not active vars
            improvement = improvement(isfinite(improvement));
            r = [r; sum(improvement)]; %#ok
            if i > 3 && sum(sign(r(end-2:end))) > 0
                solver.relaxationParameter = max(solver.relaxationParameter - 0.1, .5);
            else
                solver.relaxationParameter = min(solver.relaxationParameter + 0.1, 1);
            end
        end
        r_prev = stepReport.Residuals;
    end
    % If we converged, the last step did not solve anything
    its = i - converged;
    reports = reports(~cellfun(@isempty, reports));
    solver.relaxationParameter = omega0;
end
