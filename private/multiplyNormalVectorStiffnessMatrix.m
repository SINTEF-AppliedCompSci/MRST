function [nC,a,subGradInd] = multiplyNormalVectorStiffnessMatrix(g, C, cno, fno, nno, subhfno,weighting)
% Face-wise multiplication of stiffness matrix with normal vectors to get
% forces on the faces.
%
% This method is modified from a file in The MATLAB Reservoir Simulation
% Toolbox (MRST), see:
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

if nargin < 7 || numel(weighting) == 0
    weighting = 0;
end

[~,cn] = gridCellNodes(g,1:g.cells.num);
cn = diff(cn);
cvol = g.cells.volumes ./ cn;
% Each cell has Nd faces meeting in a vertex, meaning there will be Nd
% repetitions of [nno cno] pairs
nvol = accumarray(nno,cvol(cno)) / g.griddim;
clear cn

[a, blocksz] = rlencode([cno,nno]);
dims         = size(g.nodes.coords, 2);
assert(all(blocksz==dims));

% Expand normals, signs and permeabilities to
%  block-diagonal matrices such that nT may be constructed by matrix-matrix
%  multiplications:
[in,jn] = ndgrid(subhfno, 1:dims);
[ic,jc] = ndgrid(subhfno, 1:dims^2);
jn     = bsxfun(@plus, jn, rldecode(cumsum(blocksz)-blocksz(1), blocksz));
jc     = bsxfun(@plus, jc, rldecode(cumsum(blocksz.^2)-blocksz(1)^2, blocksz));

% Use fraction of face normal as subface normal
numnodes = double(diff(g.faces.nodePos));
N     = g.faces.normals  (fno,:)./ numnodes(fno, ones(1,dims));

N     = sparse(in,jn,N);

clear in jn

C = cat(1,C{:});
nC = cell(dims,1);
cind = (1:dims^2:max(cno) * dims^2)';
cind = reshape(bsxfun(@plus,cind,0:(dims - 1))',[],1);
% Loop over dimensions to compute the forces
for iter1 = 1 : dims
    % Pick out parts of Hook's law associated with this dimension
    Cdim = C(cind + (iter1 - 1) * dims,:);
    % Distribute cell stiffnesses onto subcells
    v = reshape(bsxfun(@minus,dims * a(:,1),fliplr(0:(dims-1)))',[],1);
    Cdim=sparse(ic,jc,reshape(Cdim(v,:)', dims^2,[])');
    if weighting
        for iter2 = 1 : dims
            Cdim(iter2:dims:end,:) = sparse(1:numel(a(:,2)),a(:,2),1) * ...
                sparse(a(:,2),1:numel(a(:,2)),cvol(a(:,1))./nvol(a(:,2))) * Cdim(iter2:dims:end,:);
        end
        
    end
    nC{iter1} = N * Cdim;
end

subGradInd = jc(1:dims:end,:);
clear ic jc C Cdim N