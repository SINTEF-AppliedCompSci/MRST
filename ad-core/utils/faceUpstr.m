function [xu, xc] = faceUpstr(flag, x, N, sz)
% Perform single-point upwinding of cell values to face
%
% SYNOPSIS:
%   [xu, xc] = faceUpstr(flag, x, N, sz)
%
% DESCRIPTION:
%   Perform single-point upwind. A robust discretization for
%   transported/hyperbolic variables.
%
% PARAMETERS:
%   flag - Boolean for each face indicating if the face should take the
%          value from the first cell (if true), or the second cell (if
%          false).
%   x    - Vector of values to be upwinded. One value per cell in the
%          domain. See `sz` input.
%   N    - Neighborship. A number of faces by 2 array. Each row corresponds
%          to the cells connected to the face with that number.
%   sz   - Vector of length 2. First entry corresponds to the number
%          of faces and the second is the total number of cells.

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


    if numel(flag) == 1
        flag = repmat(flag, size(N, 1), 1);
    end
    assert(numel(flag) == size(N, 1) && islogical(flag), ...
        'One logical upstream flag must'' be supplied per interface.');
    upCell       = N(:,2);
    upCell(flag) = N(flag,1);
    if isnumeric(x)
        % x is a simple matrix, we just extract values using the cells
        xu = x(upCell, :);
    else
        % x is likely AD, construct a matrix to achieve upstream weighting
        xu = sparse((1:sz(1))', upCell, 1, sz(1), sz(2))*x;
    end
    if nargout > 1
        xc = x;
    end
end