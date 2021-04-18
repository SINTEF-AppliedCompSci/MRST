function g = get2DGridFromFaces(G, f)
% Make a surface grid from a subset of grid faces. Since we may have
% crossing faces (>2 faces for an edge), add a sign indicator for 
% g.cells.faces. Keep same node-orientation of (2D grid) cells as the 
% original (3D grid) faces 

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

if islogical(f)
    f = find(f); 
end
[nix1, nix2, faceNo] = getFaceSegments(G, f);
% make local node-index
mapN = sparse(G.nodes.num,1);
mapN([nix1;nix2]) = 1;
nNodes = nnz(mapN);
mapN(mapN>0) = (1:nNodes)';
segs = full([mapN(nix1), mapN(nix2)]);
[~, nn] = rlencode(faceNo);

% make local cellNo (3D faceNo)
% mapF        = sparse(G.faces.num,1);
% mapF(faceNo)= 1;
% nCells     =  nnz(mapF);
% mapF(mapF>0) = (1:nCells);


cellNo = rldecode((1:numel(f))', nn);
% cut from makePlanarGrid.m --------------------------------------------- 
% Identify unique edges
[e,i]      = sort(segs, 2);
[~,j]      = sortrows(e);
k(j)       = 1:numel(j);
[edges, n] = rlencode(e(j,:));
if any(n>2)
    warning('get2DGridFromFaces needs updatating, currently assumes 2<= faces per segment ');
end
edgenum    = rldecode((1:size(edges,1))', n);
edgenum    = edgenum(k);
edgesign   = i(:,1);

% Identify neigbors assuming two per edge...
p          = sub2ind(size(edges), edgenum, edgesign);
neigh      = zeros(size(edges));
neigh(p)   = cellNo;
%------------------------------------------------------------------------
nf = size(edges,1);
g.nodes = struct('num', nf, 'coords', G.nodes.coords(mapN>0, :));
g.faces = struct('num', nf, 'nodePos', (1:2:(2*nf+1))', 'neighbors', neigh, ...
                 'tag', zeros(nf, 1), 'nodes', reshape(edges', [], 1));
g.cells = struct('num', numel(f), 'facePos', cumsum([1; nn]), 'indexMap', nan(numel(f),1), ...
                 'parent', f, 'faces', [edgenum, zeros(numel(edgenum),1)]);
g.griddim = 2;
g.type = {'sliceGrid'};
g.faces.isX = n>2;
end    

% -------------------------------------------------------------------------
function [nodes1, nodes2, faceNo] = getFaceSegments(G, f)
[np1, np2] = deal(G.faces.nodePos(f), G.faces.nodePos(f+1));
nNodes = np2-np1;
faceNo = rldecode(f(:), nNodes);
nodes1 = G.faces.nodes(mcolon(np1, np2-1));
locpos = [1; cumsum(nNodes)+1];
next   = (2:numel(nodes1)+1) .';
next(locpos(2 : end) - 1) = locpos(1 : end-1);
nodes2 = nodes1(next);
end
