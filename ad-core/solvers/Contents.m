% SOLVERS
%
% Files
%   AMGCL_CPRSolverAD      - Linear solver that calls external compiled multigrid solver
%   AMGCL_CPRSolverBlockAD - Linear solver that calls external compiled multigrid solver
%   AGMGSolverAD           - Linear solver that calls external compiled multigrid solver
%   AMGCLSolverAD          - Linear solver that calls external compiled multigrid solver
%   AMGCLSolverBlockAD     - Linear solver that calls external compiled multigrid solver
%   BackslashSolverAD      - Linear solver that calls standard MATLAB direct solver mldivide "\"
%   CPRSolverAD            - Solve a problem with a pressure component using constrained a pressure residual method
%   getNonLinearSolver     - Set up reasonable defaults for the nonlinear solver for a field
%   GMRES_ILUSolverAD      - Preconditioned GMRES solver.
%   HandleLinearSolverAD   - Simple solver for wrapping functions on the form x = fn(A, b);
%   LinearizedProblem      - A linearized problem within a non-linear iteration
%   LinearSolverAD         - Base class for linear solvers in the AD framework
%   NonLinearSolver        - Generalized Newton-like nonlinear solver
%   NoOpSolverAD           - Linear solver that does nothing.

%{
Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.

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
