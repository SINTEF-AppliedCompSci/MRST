% SOLVERS
%
% Files
%   AGMGSolverAD       - Linear solver that calls external compiled multigrid solver
%   BackslashSolverAD  - Linear solver that calls standard MATLAB direct solver mldivide "\"
%   CPRSolverAD        - Solve a problem with a pressure component using constrained a pressure residual method
%   getNonLinearSolver - Set up reasonable defaults for the nonlinear solver for a field
%   GMRES_ILUSolverAD  - Preconditioned GMRES solver.
%   LinearizedProblem  - A linearized problem within a non-linear iteration
%   LinearSolverAD     - Base class for linear solvers in the AD framework
%   NonLinearSolver    - Generalized Newton-like nonlinear solver
%   NoOpSolverAD       - Linear solver that does nothing.
%   solveMinistep      - Attempt to solve a single mini timestep while trying to avoid

%{
#COPYRIGHT#
%}
