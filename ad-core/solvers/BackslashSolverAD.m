classdef BackslashSolverAD < LinearSolverAD
    % Linear solver that calls standard MATLAB direct solver mldivide "\"
    %
    % SYNOPSIS:
    %   solver = BackslashSolverAD()
    %
    % DESCRIPTION:
    %   This solver solves linearized problems using matlab builtin mldivide.
    %
    % SEE ALSO:
    %   `LinearSolverAD`

   methods
       function solver = BackslashSolverAD(varargin)
           solver@LinearSolverAD();
           solver = merge_options(solver, varargin{:});
       end
       
       function [result, report] = solveLinearSystem(solver, A, b, varargin)
           timer = tic();
           result = A\b;
           t_solve = toc(timer);
           report = solver.getSolveReport('LinearSolutionTime', t_solve);
       end
       
        function [d, sn] = getDescription(solver)
            sn = 'mldivide';
            sn = [sn, solver.id];
            d = 'Matlab \ operator (mldivide direct solver)';
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
