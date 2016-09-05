function distFaceCell = computeDistFaceCell(g,cno,fno,nno,subhfno,eta)
% Compute distance from cell centers to face centers.
%
% This method is modified from a file in The MATLAB Reservoir Simulation
% Toolbox (MRST), see
%   mrst/modules/mpfa/computeMultiPointTrans.m
% 
%{
Parital copyright 2009-2016 SINTEF ICT, Applied Mathematics.
Partial copyright 2016, University of Bergen.

This file is part of FVBiot.

FVBiot is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

FVBiot is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
You should have received a copy of the GNU General Public License
along with this file.  If not, see <http://www.gnu.org/licenses/>.
%}

[~, blocksz] = rlencode([cno,nno]);
dims         = size(g.nodes.coords, 2);
assert(all(blocksz==dims));
%% Expand centroid differences, normals, signs and permeabilities to
%  block-diagonal matrices such that B may be constructed by matrix-matrix
%  multiplications:
[i,j] = ndgrid(subhfno, 1:dims);
j     = bsxfun(@plus, j, rldecode(cumsum(blocksz)-blocksz(1), blocksz));

etaVec = eta * ones(numel(fno),1);
isbnd = any(g.faces.neighbors == 0,2);
etaVec(isbnd(fno)) = 0;

cp = g.faces.centroids(fno,:) + bsxfun(@times,etaVec,(g.nodes.coords(nno,:) - g.faces.centroids(fno,:)));

distFaceCell     = cp - g.cells.centroids(cno,:);

distFaceCell = sparse(i,j,distFaceCell);
clear i j cp etaVec isbnd