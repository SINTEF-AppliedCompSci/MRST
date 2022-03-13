function [A, rhs] = getSystemIncompTPFA(state, G, T, fluid, varargin)
% Extract linear system from two point solver without solving any linear
% systems.

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
    s = incompTPFA(state, G, T, fluid, ...
                                      'use_trans', numel(T) == G.faces.num, ...
                                      varargin{:}, ...
                                      'LinSolve', @(A, x) zeros(size(x, 1), 1), 'MatrixOutput', true);
    A = s.A;
    rhs = s.rhs;
end
