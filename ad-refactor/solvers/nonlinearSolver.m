classdef nonlinearSolver
    properties
        maxIterations
        maxSubsteps
        linearSolver
        verbose
    end
    
    methods
        function solver = nonlinearSolver(varargin)
            opt = struct('maxIterations', 25, ...
                         'maxSubsteps',   32, ...
                         'verbose',       mrstVerbose(), ...
                         'isAdjoint',     false, ...
                         'linearSolver',   []);
            opt = merge_options(opt, varargin{:});
            
            solver.maxIterations = opt.maxIterations;
            solver.maxSubsteps   = opt.maxSubsteps;
            solver.verbose     = opt.verbose;
            
            if isempty(opt.linearSolver)
                solver.linearSolver = mldivideSolverAD();
            else
                solver.linearSolver = opt.linearSolver;
            end
        end
        
        function [state, status] = solveTimestep(solver, state0, dT, model, varargin)
            % Solve a timestep for a non-linear system using one or more substeps
            drivingForces = struct('Wells', [],...
                                   'bc',    [],...
                                   'src',   []);
            
            drivingForces = merge_options(drivingForces, varargin{:});
            
            
            converged = false;
            done = false;
            
            dt = dT;
            ministepNo = 1;
            
            while ~done
                state = state0;
                for iter = 1:ministepNo
                    % Do a bunch of ministeps
                    state0_local = state;
                    for i = 1:solver.maxIterations
                        [state, converged] = model.stepFunction(state, state0_local, dt, drivingForces, solver.linearSolver, 'iteration', i);
                        if converged
                            break
                        end
                    end
                    if ~converged
                        break
                    end
                end
                
                if converged
                    done = true;
                elseif ministepNo < solver.maxSubsteps
                    % We didn't converge, but we are still away from the
                    % maximum number of substeps. Double the amount for the
                    % next iteration.
                    ministepNo = 2*ministepNo;
                    dt = dt/2;
                else
                    warning('Solver did not converge')
                    break
                end
            end
            
            status = struct('iterations',   iter,...
                            'converged',    converged,...
                            'ministeps', 	ministepNo);

        end
    end
end