function [state, converged, failure, its, reports] = solveMinistep(solver, model, state, state0, dt, drivingForces)
% Attempt to solve a single mini timestep while trying to avoid
% stagnation or oscillating residuals.
%
% Note: This is an internal function for the NonLinearSolver class that is
% exposed for testing purposes. Do not use directly, as the interface may
% disappear or change without prior notice. Instead, use
% NonLinearSolver.solveTimestep() instead.

%{
Copyright 2009-2016 SINTEF ICT, Applied Mathematics.

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
                dispif(solver.verbose > 0, ...
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