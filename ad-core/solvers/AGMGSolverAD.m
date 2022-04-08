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
       singleApply = false; % Just apply the preconditioner once, do not solve to tolerance
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
       
       function [result, report] = solveLinearSystem(solver, A, b, varargin)
           cleanAfter = false;
           if ~solver.setupDone
               solver.setupSolver(A, b);
               cleanAfter = true;
           end
           % Solve the linear system to a given tolerance
           if solver.reuseSetup
               code = 2 + double(solver.singleApply);
           else
               assert(~solver.singleApply, 'Cannot use singleApply without reuseSetup.');
               code = [];
           end
           [result, flag, relres, iter, resvec] = agmg(A, b,...
                                                       [],... % Restart option (default)
                                                       solver.tolerance, ...
                                                       solver.maxIterations, ...
                                                       double(solver.verbose > 1), ... % If verbose > 1, do internal verbose
                                                       [], ...
                                                       code);
           report = solver.getSolveReport('Converged', flag < 1, ...
                                          'Residual', relres, ...
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
               if solver.setupDone
                   solver = solver.cleanupSolver(A, b);
               end
               dispif(solver.verbose, 'Setting up AGMG preconditioner...');
               v = double(solver.verbose > 1);
               agmg(A, [], [], [], [], v, [], 1);
               solver.setupDone = true;
               dispif(solver.verbose, ' AGMG set up.\n');
           end
       end
       
       function  solver = cleanupSolver(solver, A, b, varargin)
           % Clean up solver after use (if needed)
           if solver.reuseSetup && solver.setupDone
               dispif(solver.verbose, 'Cleaning up AGMG...');
               v = double(solver.verbose > 1);
               agmg(A, [], [], [], [], v, [], -1);
               solver.setupDone = false;
               dispif(solver.verbose, ' AGMG reset.\n');
           end
       end

       function [d, sn] = getDescription(solver)
           sn = 'AGMG';
           sn = [sn, solver.id];
           d = 'Aggregation based multigrid (AGMG)';
       end
       
       function delete(solver)
           solver.cleanupSolver(solver, [], []);
       end
   end
end

%{
Copyright 2009-2022 SINTEF Digital, Mathematics & Cybernetics.

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
