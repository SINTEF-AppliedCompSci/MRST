function iA = invertDiagonalBlocks(A, sz)
% Invert block diagonal matrices by Matlab.
%
% This method is identical to a subroutine in MRST, 
%   see mrst/modules/mpfa/computeMultiPointTrans.m
% It has been copied to a separate file to enable access from other
% funcions.
%
% Copyright statement from the original file is pasted below.
%{
Copyright 2009-2016 SINTEF ICT, Applied Mathematics.

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

V = zeros([sum(sz .^ 2), 1]);
[p1, p2] = deal(0);

for b = 1 : numel(sz),
    n  = sz(b);
    n2 = n * n;
    i  = p1 + (1 : n);
    
    V(p2 + (1 : n2)) = inv(full(A(i, i)));
    
    p1 = p1 + n;
    p2 = p2 + n2;
end

iA = blockDiagMatrix(V, sz);
clear V;