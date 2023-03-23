function [cutFaces, cutCells, coords, sliceFaces] = extractGridInfoFromPolygons(G, p)
% Identyfy unique coords, later remove/ignore those coinciding with grid nodes

%{
Copyright 2009-2023 SINTEF Digital, Mathematics & Cybernetics.

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

[coords, ~, uix] = uniquetol(p.coords3D, sqrt(eps), 'ByRows', true);
uix   = uix'; % case single poly
p.nodes = uix(:);
% rotate nodes such that equal node numbers are adjacent (e.g., not first and last)
p = rotatePolygonNodes(p);
sliceFaces = [];

% cut-cells are those beeing split by non-degenerate polygons
nPoly = numel(p.nodePos)-1; 
if nPoly == 0
     coords = [];
     [cutFaces, cutCells] = deal(struct('ix', []));
    return
end
% index to first appearence of node
uniqueIx = [true; diff(p.cellIx) | diff(p.nodes)]; 
polyNo   = rldecode( (1:nPoly)', diff(p.nodePos));
nNodes   = accumarray(polyNo, uniqueIx);
% disregard degenerates
validPoly  = nNodes > G.griddim -1;
uniqueIx   = uniqueIx & validPoly(polyNo);
% cellIx per polygon
cellIx   = p.cellIx(p.nodePos(validPoly));   
% update
nNodes  = nNodes(validPoly);
% unique represenation for polygon inside cut-cell
nodes_cut  = p.nodes(uniqueIx);
% compute additional required geometry (not needed in case of topo split)
if G.griddim == 3
    [normals, areas] = computeCutCellNormals(coords, nodes_cut, nNodes);
elseif G.griddim == 2
    % Assume for now planar grid. Need normal to plane
    pn = getPlaneNormal(G);
    %dispif(mrstVerbose, 'WARNING: assuming 2D-grid is planar, needs updating\n%s', ...
    %                    '(extractGridInfoFromPlygons.m)');
    [normals, areas] = computeCutCellNormals(coords, nodes_cut, nNodes, pn);
end

cutCells = struct('ix',       cellIx, ...
                  'nNodes',   nNodes, ...
                  'nodes',    nodes_cut, ...
                  'normals',  normals, ...
                  'areas',    areas);

%--------------------------------------------------------------------------
% cut-faces are those beeing split by line-segments not along face
% segments
% find face-node pairs
[sr, ia] = unique([p.faceIx, p.nodes], 'rows');
[faceIx, nodes] = deal(sr(:,1), sr(:,2));      
% update
localNodePos = p.locPos(ia);
tValue       = p.tValue(ia);
% OK with strict tollerance here (end t-values should be exactly 0/1)
onNodeTol  = sqrt(eps);
onGridNode = tValue < onNodeTol | tValue > 1 - onNodeTol;

% collect unique faces and number of occurances (2 for 3D, 1 for 2D)
[fix, nn] = rlencode(faceIx);
fixPos    = cumsum([1;nn]);

% keep track of which faces needs updating, and coords coincinding with 
% grid-nodes. A face could need updating without beeing split, e.g., just
% add a single node.
isUpdateFace    = false(size(fix));
isSplitFace     = false(size(fix));

undefinedCount = 0;
coordOnNode     = nan(size(coords,1), 1);
coordOnNode(nodes) = p.snapNode(ia); 
for k = 1:numel(fix)
    if nn(k) == G.griddim-1
        curIx = fixPos(k):fixPos(k+1)-1;
        if all(~onGridNode(curIx))
            isSplitFace(k) = true;
        else
            % if a valid split, there should be at least one face-node
            % between the two intersection points (the oppposite can happen
            % as a result of snapping)
            tv = tValue(curIx);
            onNode  = onGridNode(curIx);
            if any(~onNode)
                isUpdateFace(k) = true;
            end
            nodePos = G.faces.nodePos(fix(k) + [0 1]');
            fnodes  = G.faces.nodes(nodePos(1):(nodePos(2)-1));
            if G.griddim == 3
                nf     = diff(nodePos);
                % intersection is between node np and next face-node with
                % relative distance tv
                % Give each intersection point a position relative to local
                % face node number, such that position 1 is between first
                % and last, 2 is on first, 3 is between first and second
                % etc.
                expr = round(2*localNodePos(curIx) + ~onNode + 2*tv.*onNode);
                locpos = mod(expr-1, 2*nf) + 1;
                % to be a valid cut, there needs to be at least one node
                % between the two intersections, hence distance in locpos
                % greater than to 2 (since at least one of the intersection 
                % points lie on a node)
                dist   = min(mod([diff(locpos), diff(locpos([2,1]))], 2*nf));
                if dist > 2
                    isSplitFace(k) = true;
                end
                % identify coinciding node(s)
                for kc = 1:2
                    if onNode(kc)
                        coordno = nodes(curIx(kc));
                        coordOnNode(coordno) = fnodes(locpos(kc)/2);
                    end
                end
            elseif G.griddim == 2
                if round(tv) == 0
                    coordno = nodes(curIx);
                    coordOnNode(coordno) = fnodes(1);
                end
            end
        end
    elseif nn(k) >= 2
        undefinedCount = undefinedCount +1;
    end
end
if undefinedCount > 0
    dispif(mrstVerbose, 'Found %d undefined face-splits, assumed co-planar.\n', ...
            undefinedCount);
end

% if there are any coords on grid-nodes, we need to remove the corresponding
% poly-nodes and re-index poly-nodes
nodeReIndex = (1:size(coordOnNode))';
if ~all(isnan(coordOnNode))
    coix   = ~isnan(coordOnNode);
    nodeNo = coordOnNode(coix);
    % we could remove these coords, but don't bother now (stowaways)
    % use convention that negative node-no points to grid node, positive to
    % (new) poly-node
    nodeReIndex(coix) = -nodeNo;
end
    
% re-index for cut-cells
cutCells.nodes = nodeReIndex(cutCells.nodes);

% make cut-faces
ix         = cumsum([1;nn]);
updIx      = isSplitFace | isUpdateFace;
nonSplitIx = ~isSplitFace & isUpdateFace;
ix         = ix(updIx);
if G.griddim == 3
    cutFaces = struct('ix',          fix(updIx), ...
                      'nodes',       nodeReIndex([nodes(ix), nodes(ix+1)]), ...
                      'nodePos',     [localNodePos(ix), localNodePos(ix+1)], ...
                      'isTrueSplit', ~nonSplitIx(updIx));
elseif G.griddim == 2
    cutFaces = struct('ix',           fix(isSplitFace), ...
                      'nodes',        nodeReIndex(nodes(ix)), ...
                      'nodePos',      localNodePos(ix), ...
                      'isTrueSplit', ~nonSplitIx(updIx));
end
end

%----------------------------------------------------------------------
function [n, area] = computeCutCellNormals(coords, nodes, nNodes, normal)
griddim = 3 - (nargin == 4);
nc = numel(nNodes);
if nc > 0
    if griddim == 3
        cNo   = rldecode((1:nc)', nNodes);
        acc   = sparse(cNo, (1:numel(cNo))', 1);
        cent  = (acc*coords(nodes,:))./nNodes;
        %sub-triangle normals
        pos  = cumsum([1;nNodes]);
        next = getNext(nodes, pos);
        [a,b] = deal(coords(nodes,:), coords(next,:));
        nt = cross(b-a, cent(cNo,:)-a, 2)./2;
        na = (acc*nt);
    else % for 2D, use plane-normal, this must be updated for curved grids!
        normal = normal(:)/norm(normal);
        [a,b] = deal(coords(nodes(1:2:end-1),:), coords(nodes(2:2:end),:));
        na = cross(b-a, ones(numel(nodes)/2,1)*normal', 2);
    end
    area = sqrt(sum(na.^ 2, 2));
    n = bsxfun(@rdivide, na, area);
else
    [n, area] = deal([]);
end
end

%----------------------------------------------------------------------
function pn = getPlaneNormal(G)
    if isfield(G.cells, 'normal')
        pn = G.cells.normal(1,:)';
    else
        % find normal of any cell
        if ~isfield(G.cells, 'volumes')
            G = computeGeometry(G);
        end
        [~, c] = max(G.cells.volumes);
        fNo = G.cells.facePos(c):(G.cells.facePos(c+1)-1);
        f   = G.cells.faces(fNo);
        v   = diff(G.faces.centroids(f,:));
        pn  = sum(cross(v(1:end-1,:), v(2:end,:),2));
        pn = pn/norm(pn);
    end
end
        
%----------------------------------------------------------------------
function next = getNext(n, p)
next               = (2 : size(n, 1) + 1) .';
next(p(2:end) - 1) = p(1 : end-1);
next = n(next);
end
%----------------------------------------------------------------------
%{
function ncan = canonicalPoly(n)
% remove repeated and rotate to smallest index first
np  = size(n,1);
c   = cell(np,1);
uix = [true(np, 1), abs(diff(n, [], 2))>0];
maxn = max(sum(uix,2));
for k = 1:np
    p = n(k, uix(k,:));
    [~, ix] = min(p);
    p = p([ix:end, 1:ix-1]);
    p = [p, p(end)*ones(1, maxn-numel(p))]; %#ok (not the case)
    c{k} = p;
end
ncan = vertcat(c{:});
end
%}
%----------------------------------------------------------------------
function p = rotatePolygonNodes(p)
    % check if first and last node are equal
    ii = find( p.nodes(p.nodePos( 1:(end-1) )) == ...
               p.nodes(p.nodePos(2:end)-1) );
    if ~isempty(ii)
        ni = numel(ii);
        [ix, rix] = deal(cell(ni, 1));
        for k = 1:ni
            ix{k} = (p.nodePos(ii(k)) : p.nodePos(ii(k)+1)-1)';
            n1 = p.nodes(ix{k}(1));
            for shift = numel(ix{k})-1:-1:1
                if p.nodes(ix{k}(shift))~=n1
                    break
                end
            end
            if ~(shift==1)
                rix{k} = ix{k}([shift+1:end, 1:shift]);
            else % all nodes equal
                [rix{k}, ix{k}] = deal([]);
            end
        end
        [ix, rix] = deal(vertcat(ix{:}), vertcat(rix{:}));
        % rotate index of all relevant fields
        fn = fieldnames(p);
        for k = 1:numel(fn)
            if size(p.(fn{k}),1) == numel(p.nodes)
                p.(fn{k})(ix) = p.(fn{k})(rix);
            end
        end
    end
end

%{
WORK IN PROGRESS!
function p = processNonConformalSnapping(G, p, uniqueIx, nNodes, cellIx)
% To be valid, a cut can maximally share one segment with each cell-face or
% be equal to one of the cell faces.
% Only consider those with >= 3 snapped nodes
polyNo = rldecode( (1:numel(nNodes))', nNodes);
nSnap    = accumarray(polyNo, ~isnan(p.snapNode(uniqueIx)));
checkIx  = find(nSnap >= 3);
if ~isempty(checkIx)
    facePos = G.cells.facePos;
    nodePos = G.faces.nodePos;
    polyNodes = something;
    for k = 1:numel(checkIx)
        c  = cellIx(checkIx(k));
        fp = facePos(c):(facePos(c+1)-1);
        f  = G.cells.faces(fp);
        for kf = 1:numel(f)
            np = nodePos(f(kf)):(nodePos(f(kf)+1)-1);
            n  = G.faces.nodes(np);
            isSame = ismember(n, polyNodes);
            if nnz(isSame) >= 3
                i1 = 1;
                for kn = 1:numel(n)
                    if isSame(kn)
                        i2 = kn;
                    else
                        if i2-i1 >= 2
                            % we have >= three consecutive
                            % remove intermediate nodes from cutcell
                            i1+1:i2-1
                            % add intermediate nodes to cutFaces
                            i1:i2 
                        else
                            [i1, i2] = deal(kn);
                        end
                    end
                end
            end
        end
    end
end
end
%}                            
