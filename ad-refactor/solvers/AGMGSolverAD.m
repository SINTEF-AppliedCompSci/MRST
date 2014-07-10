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
%   BackslashSolverAD

%{
Copyright 2009-2014 SINTEF ICT, Applied Mathematics.

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
       
   end
   methods
       function solver = AGMGSolverAD(varargin)
            require agmg
            solver = solver@LinearSolverAD(varargin{:});
       end
       
       function [result, report] = solveLinearSystem(solver, A, b)
           % Solve the linear system to a given tolerance
           [result, flag, relres, iter, resvec] = agmg(A, b, [],...
                        solver.tolerance, solver.maxIterations, [], [], 2);
           report = struct('Converged', flag < 1, ...
                           'RelativeResidual', relres, ...
                           'Iterations',   iter);
            if solver.extraReport
                report.ResidualHistory = resvec;
            end
       end
       
       function solver = setupSolver(solver, A, b, varargin) %#ok 
           % Run setup on a solver for a given system
           agmg(A,[],[],[],[],[],[], 1);
       end
       
       function  solver = cleanupSolver(solver, A, b, varargin) %#ok 
           % Clean up solver after use (if needed)
           agmg(A,[],[],[],[],[],[], -1);
       end
   end
end
