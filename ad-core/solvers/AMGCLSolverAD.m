classdef AMGCLSolverAD < LinearSolverAD
% Linear solver that calls external compiled multigrid solver
%
% SYNOPSIS:
%   solver = AMGCLSolverAD()
%
% DESCRIPTION:
%    AD-interface for the AMGCL interface.
%
% NOTE:
%    This solver requires AMGCL to be installed and working.
%
% SEE ALSO:
%   BackslashSolverAD

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
   properties
       coarsening
       solver
       relaxation
       reduceToCell
       preconditioner
   end
   methods
       function solver = AMGCLSolverAD(varargin)
            require linearsolvers
            solver = solver@LinearSolverAD();
            solver.preconditioner = 'amg';
            solver.coarsening = 'smoothed_aggregation';
            solver.relaxation = 'spai0';
            solver.solver     = 'bicgstab';
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
            report = struct();
            result = callAMGCL(A, b, ...
                 'coarsening',     solver.coarsening, ...
                 'preconditioner', solver.preconditioner, ...
                 'relaxation',     solver.relaxation, ....
                 'solver',         solver.solver,...
                 'maxIterations',  solver.maxIterations, ...
                 'tolerance',      solver.tolerance);
       end
   end
end
