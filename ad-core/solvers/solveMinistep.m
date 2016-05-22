function [state, converged, failure, its, reports] = solveMinistep(solver, model, state, state0, dt, drivingForces)
% Attempt to solve a single mini timestep while trying to avoid
% stagnation or oscillating residuals.
   omega0 = solver.relaxationParameter;
   if model.stepFunctionIsLinear
       maxIts = 0;
   else
       maxIts = solver.maxIterations;
   end
   reports = cell(maxIts, 1);
   for i = 1:(maxIts + 1)
      % If we are past maximum number of iterations, step function will
      % just check convergence and return
      [state, stepReport] = ...
          model.stepFunction(state, state0, dt, drivingForces, ...
                             solver.LinearSolver, solver, ...
                             i);
      converged  = stepReport.Converged;
      failure    = stepReport.Failure;
      reports{i} = stepReport;
      if converged
         break
      end
      if failure
         break
      end

      if i > 1 && solver.enforceResidualDecrease
         if all(stepReport.Residuals >= prev_best)
            % We are not seeing reduction, but rather increase in the
            % residuals. Break and let the solver decide to either
            % abort or cut the timestep.
            break;
         end
      end
      prev_best = stepReport.Residuals;

      if solver.useRelaxation
         % Store residual history during nonlinear loop to detect
         % stagnation or oscillations in residuals.
         if i == 1
            res = nan(maxIts + 1, numel(stepReport.Residuals));
         end
         res(i, :) = stepReport.Residuals;

         isOk = res(i, :) <= model.nonlinearTolerance;
         isOscillating = solver.checkForOscillations(res, i);
         isStagnated = solver.checkForStagnation(res, i);
         % We will use relaxations if all non-converged residuals are
         % either stagnating or oscillating.
         bad = (isOscillating | isStagnated) | isOk;
         relax = all(bad) && ~all(isOk);
         if relax
            dispif(solver.verbose, ...
                   'Convergence issues detected. Applying relaxation to Newton solver.\n');
            solver.relaxationParameter = max(solver.relaxationParameter - solver.relaxationIncrement, solver.minRelaxation);
         else
            solver.relaxationParameter = min(solver.relaxationParameter + solver.relaxationIncrement, solver.maxRelaxation);
         end
      end
   end
   % If we converged, the last step did not solve anything
   its = i - converged;
   reports = reports(~cellfun(@isempty, reports));
   solver.relaxationParameter = omega0;
   if converged
      [state, r] = model.updateAfterConvergence(state0, state, dt, drivingForces);
      reports{end}.FinalUpdate = r;
   end
end