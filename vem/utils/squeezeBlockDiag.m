function A = squeezeBlockDiag(A, n, r, c)
%Squeezes a block diagonal matrix in which each block has the same number
%of columns OR the same numer of rows.
%
%   The result is a 'column' or 'row' matrix:
%
%   [A     ]             [A]         [A     ]
%   [   B  ] -{column}-> [B]   and   [   B  ] -{row}-> [A  B  C]
%   [     C]             [C]         [     C]
%
%   SYNOPSIS:
%       A = squeezeBlockDiag(A, n, r, c)
%
%   REQUIRED PARAMETERS:
%       A - Block diagonal matrix in which each block has the same number
%           of columns OR the same numer of rows.
%
%       n - Number of matrix blocks.
%
%       r - Number of rows in resulting matrix.
%
%       c - Number of columns in resulting matrix.
%
%   RETURNS:
%       A - Squeezed matrix.

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

%   Written by Øystein Strengehagen Klemetsdal, SINTEF/NTNU, 2016.

[R, ~] = size(A);
N = numel(n);

if r == R
    A = A*repmat(speye(c),N,1);
else
    A = repmat(speye(r),1,N)*A;
end

end
