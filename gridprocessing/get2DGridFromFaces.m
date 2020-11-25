function g = get2DGridFromFaces(G, f)
if islogical(f)
    f = find(f); 
end
[nix1, nix2, faceNo] = getFaceSegments(G, f);
% make local node-index
map = sparse(G.nodes.num,1);
map([nix1;nix2]) = 1;
nNodes = sum(map);
map(map>0) = (1:nNodes)';
segs = [map(nix1), map(nix2)];
[~, nn] = rlencode(faceNo);

cellNo = rldecode((1:numel(f))', nn);
% cut from makePlanarGrid.m --------------------------------------------- 
% Identify unique edges
[e,i]      = sort(segs, 2);
[~,j]      = sortrows(e);
k(j)       = 1:numel(j);
[edges, n] = rlencode(e(j,:));
edgenum    = rldecode((1:size(edges,1))', n);
edgenum    = edgenum(k);
edgesign   = i(:,1);

% Identify neigbors assuming two per edge...
p          = sub2ind(size(edges), edgenum, edgesign);
neigh      = zeros(size(edges));
neigh(p)   = cellNo;
%------------------------------------------------------------------------
nf = size(edges,1);
g.nodes = struct('num', nf, 'coords', G.nodes.coords(map>0, :));
g.faces = struct('num', nf, 'nodePos', (1:2:(2*nf+1))', 'neighbors', neigh, ...
                 'tag', zeros(nf, 1), 'nodes', reshape(edges', [], 1));
g.cells = struct('num', numel(f), 'facePos', cumsum([1; nn]), 'indexMap', nan(numel(f),1), ...
                 'parent', f, 'faces', [edgenum, zeros(numel(edgenum),1)]);
g.griddim = 2;
g.type = {'sliceGrid'};
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