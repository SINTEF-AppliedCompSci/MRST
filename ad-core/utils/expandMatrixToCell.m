function c = expandMatrixToCell(matrix, subset)
% Expand a matrix into cell arrays. Typical usage: Converting state
% representation of composition (as matrix) into AD-values (as cell array
% of columns vectors).

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

    if iscell(matrix)
        c = matrix;
        return
    end
    if nargin == 1
        subset = ':';
    end
    
    n = size(matrix, 2);
    c = cell(1, n);
    for i = 1:n
        c{i} = matrix(subset, i);
    end
end
