classdef NoOpSolverAD < LinearSolverAD
    % Linear solver that does nothing.
    %
    % SYNOPSIS:
    %   solver = NoOpSolverAD()
    %
    % DESCRIPTION:
    %   Debug solver. It has the correct interfaces, but it always returns zero
    %   as the solution for the problem. It is, however, very fast...
    %
    % NOTE:
    %   You should not use this solver.
    %
    % SEE ALSO:
    %   BackslashSolverAD

   properties
       
   end
   methods
       
       function [result, report] = solveLinearSystem(solver, A, b)
           result = zeros(size(b));
           report = struct('Converged', false);
       end
   end
end

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
