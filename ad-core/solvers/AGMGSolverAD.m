classdef AGMGSolverAD < LinearSolverAD
    % Linear solver that calls external compiled multigrid solver
    %
    % SYNOPSIS:
    %   solver = AGMGSolverAD()
    %
    % DESCRIPTION:
    %   This solver calls the AGMG package and supports setup/cleanup of the
    %   multigrid solver for multiple solves of the same problem.
    %
    % NOTE:
    %    This solver requires AGMG to be installed and working.
    %
    % SEE ALSO:
    %   `BackslashSolverAD`

   properties
       setupDone % Internal book-keeping variable
       reuseSetup % Will reuse the setup phase to improve speed for e.g. a GMRES loop with the same matrix system. However, some systems report segfaults with this option enabled.
   end

   methods
       function solver = AGMGSolverAD(varargin)
            require agmg
            solver = solver@LinearSolverAD();
            solver.setupDone = false;
            solver.reuseSetup = false;
            solver.reduceToCell = true;
            solver = merge_options(solver, varargin{:});
       end
       
       function [dx, result, report] = solveLinearProblem(solver, problem, model)
           if solver.reduceToCell
               [dx, result, report] = solver.solveCellReducedLinearProblem(problem, model);
           else
               [dx, result, report] = solveLinearProblem@LinearSolverAD(solver, problem, model);
           end
       end
       
       function [result, report] = solveLinearSystem(solver, A, b)
           cleanAfter = false;
           if ~solver.setupDone
               solver.setupSolver(A, b);
               cleanAfter = true;
           end
           % Solve the linear system to a given tolerance
           if solver.reuseSetup
               fn = @(A, b) agmg(A, b, [], solver.tolerance, ...
                                    solver.maxIterations, [], [], 1);
           else
               fn = @(A, b) agmg(A, b, [], solver.tolerance, ...
                                    solver.maxIterations);
           end
           [result, flag, relres, iter, resvec] = fn(A, b);
           report = struct('Converged', flag < 1, ...
                           'RelativeResidual', relres, ...
                           'Iterations',   iter);
            if solver.extraReport
                report.ResidualHistory = resvec;
            end
            
            if cleanAfter
                solver.cleanupSolver(A, b);
            end
       end
       
       function solver = setupSolver(solver, A, b, varargin)
           % Run setup on a solver for a given system
           if solver.reuseSetup
               agmg(A,[],[],[],[],[],[], 1);
               solver.setupDone = true;
           end
       end
       
       function  solver = cleanupSolver(solver, A, b, varargin)
           % Clean up solver after use (if needed)
           if solver.reuseSetup
               agmg(A,[],[],[],[],[],[], -1);
               solver.setupDone = false;
           end
       end
   end
end

%{
Copyright 2009-2017 SINTEF ICT, Applied Mathematics.

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
