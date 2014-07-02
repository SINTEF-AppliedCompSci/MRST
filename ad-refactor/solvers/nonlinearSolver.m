classdef nonlinearSolver
    properties
        maxIterations
        maxSubsteps
        linearSolver
        verbose
        isAdjoint
    end
    
    methods
        function solver = nonlinearSolver(varargin)            
            solver.maxIterations = 25;
            solver.maxSubsteps   = 25;
            solver.verbose       = mrstVerbose();
            solver.maxSubsteps   = 32;
            solver.isAdjoint     = false;
            solver.linearSolver  = [];
            
            solver = merge_options(solver, varargin{:});
            
            if isempty(solver.linearSolver)
                solver.linearSolver = mldivideSolverAD();
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
            ministepNo = 1;
            
            itCount = 0;
            while ~done
                state = state0;
                
                reports = cell(ministepNo, 1);
                for iter = 1:ministepNo
                    % Do a bunch of ministeps
                    state0_local = state;
                    
                    nonlinearReports = cell(solver.maxIterations, 1);
                    for i = 1:solver.maxIterations
                        [state, stepReport] = ...
                            model.stepFunction(state, state0_local, dt, drivingForces, ...
                                               solver.linearSolver, 'iteration', i);
                        itCount = itCount + 1;
                        converged  = stepReport.Converged;

                        if converged
                            break
                        end
                        nonlinearReports{i} = stepReport;
                    end
                    if ~converged
                        break
                    end
                    reports{iter} = struct();
                    reports{iter}.NonlinearReport = nonlinearReports(~cellfun(@isempty, nonlinearReports));
                    reports{iter}.LocalTime = dt*iter;
                    % This line does nothing atm
                    reports{iter}.Converged = converged;
                end
                
                
                if converged
                    done = true;
                elseif ministepNo < solver.maxSubsteps
                    % We didn't converge, but we are still away from the
                    % maximum number of substeps. Double the amount for the
                    % next iteration.
                    warning('Solver did not converge, cutting timestep')
                    ministepNo = 2*ministepNo;
                    dt = dt/2;
                else
                    warning('Solver did not converge')
                    break
                end
            end
            % Truncate reports from step functions
            reports = reports(~cellfun(@isempty, reports));
            report = struct('Iterations',       itCount,...
                            'Converged',        converged,...
                            'MinistepCount', 	ministepNo);
            % Add seperately because struct constructor interprets cell
            % arrays as repeated structs.
            report.StepReports = reports;
        end
    end
end