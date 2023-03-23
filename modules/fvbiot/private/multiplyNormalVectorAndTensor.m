function [nK,a,j] = multiplyNormalVectorAndTensor(g, T, cno, fno, nno, subhfno)
% Multiply normal vectors and permeability tensors on a sub-cell level
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


[a, blocksz] = rlencode([cno,nno]);
dims         = size(g.nodes.coords, 2);
assert(all(blocksz==dims));

% Expand normals, signs and permeabilities to
%  block-diagonal matrices such that nT may be constructed by matrix-matrix
%  multiplications:
[i,j] = ndgrid(subhfno, 1:dims);
j     = bsxfun(@plus, j, rldecode(cumsum(blocksz)-blocksz(1), blocksz));

% Use fraction of face normal as subface normal
numnodes = double(diff(g.faces.nodePos));
N     = g.faces.normals  (fno,:)./ numnodes(fno, ones(1,dims));
N     = sparse(i,j,N);

K     = sparse(i,j,reshape(T(a(:,1),:)', dims, [])');
nK = N * K;
j = j(1:dims:end,:);