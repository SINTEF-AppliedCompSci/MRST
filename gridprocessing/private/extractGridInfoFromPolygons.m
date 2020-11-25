function [cutFaces, cutCells, coords, sliceFaces] = extractGridInfoFromPolygons(G, p)
% Identyfy unique coords, later remove/ignore those coinciding with grid nodes
[coords, ~, uix] = uniquetol(p.coords3D, sqrt(eps), 'ByRows', true);
uix   = uix'; % case single poly
% map to unique coords, keep original for indexing into tValue,faceIx etc.
nodesMUC = uix(p.nodes);
nodesM   = p.nodes;
cellIx   = p.cellIx;
%check for equal polygons modulo repeated nodes and rotation of vertices
[~, ia] = unique(canonicalPoly(nodesMUC), 'rows');
if numel(ia) ~= size(nodesM,1)
    nodesMUC  = nodesMUC(ia,:);
    nodesM    = nodesM(ia,:);
    cellIx    = p.cellIx(ia,:);
end
sliceFaces = [];
%{
Earlier version did not work robustly, forward identifying slice-faces to
cutGrid 
--------
sliceFaces = [];
cellIx     = p.cellIx;
% two equal polygons indicate a slice along a face, remove these, but
% record in sliceFaces
if numel(ia) ~= numel(ib)
     npol = accumarray(ib, ones(size(ib)));
%     if any(np>2)
%         warning('Problematic face: >2 neighbors??')
%     end
     keep = npol(ib) == 1;
%     [~, pairOrder] = sort(ib(~keep));
%     pairOrder = [pairOrder(1:2:end), pairOrder(2:2:end)];
%     cellPairs = cellIx(~keep);
%     cellPairs = cellPairs(pairOrder);
%     sliceFaces = facesFromNeighbors(G, cellPairs);
    % update poly-info
    nodesMUC  = nodesMUC(keep,:);
    nodesM   = nodesM(keep,:);
    cellIx = cellIx(keep,:);
    % update additional fields
   % flds = {'nodes'};
   % for k  =1: numel(flds)
   %     p.(flds{k}) = p.(flds{k})(keep,:);
   % end
end
%}
%--------------------------------------------------------------------------

% cut-cells are those beeing split by non-degenerate polygons
nPoly = size(nodesM, 1);
if nPoly == 0
     [coords, g2] = deal([]);
     [cutFaces, cutCells] = deal(struct('ix', []));
    return
end
% index to first appearance of node
uniqueIx = [true(nPoly, 1), abs(diff(nodesMUC, [], 2))>0];
nNodes = sum(uniqueIx, 2);
% disregard degenerates
validPoly  = nNodes > G.griddim -1;
% update
nodesMUC    = nodesMUC(validPoly, :);
nodesM     = nodesM(validPoly, :);
uniqueIx   = uniqueIx(validPoly, :);
% vectorize
nodesUC  = reshape(nodesMUC', [], 1);
cellIx   = cellIx(validPoly);   
% update
nNodes  = nNodes(validPoly);

% unique represenation for polygon inside cut-cell
nodes_cut  = nodesUC(uniqueIx');
% compute additional required geometry
if G.griddim == 3
    [normals, areas] = computeCutCellNormals(coords, nodes_cut, nNodes);
elseif G.griddim == 2
    % Assume for now planar grid. Need normal to plane
    pn = getPlaneNormal(G);
    dispif(mrstVerbose, 'WARNING: assuming 2D-grid is planar, needs updating\n%s', ...
                        '(extractGridInfoFromPlygons.m)');
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
% local face-node index corresponding to start of cut segment
locIx = p.nodeIx(:,1); 
% expand to each polygon segment
localNodePos  = reshape(locIx(nodesM)', [], 1);
tValue        = reshape(p.tValue(nodesM)', [], 1);
faceIx        = reshape(p.faceIx(nodesM)', [], 1);


% find face-node pairs
[sr, ia]    = unique([faceIx, nodesUC], 'rows');
[faceIx, nodesUC] = deal(sr(:,1), sr(:,2));         % reset nodesUC here?
% update
localNodePos = localNodePos(ia); 
tValue       = tValue(ia);

% OK with strict tollerance here (end t-values should be exactly 0/1)
onNodeTol    = sqrt(eps);
onGridNode = tValue < onNodeTol | onNodeTol > 1 - sqrt(eps);

% collect unique faces and number of occurances (2 for 3D, 1 for 2D)
[fix, nn] = rlencode(faceIx);
fixPos    = cumsum([1;nn]);

% keep track of split-faces, problems and coords coincinding with grid-nodes
isSplitFace     = false(size(fix));
unexpectedCount = 0;
coordOnNode     = nan(size(coords,1), 1);
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
                % to by a valid cut, there needs to be at least one node
                % between the two intersections, hence distance in locpos
                % greater than or equal to 2
                dist   = min(mod([diff(locpos), diff(locpos([2,1]))], 2*nf));
                if dist >= 2
                    isSplitFace(k) = true;
                end
                % identify coinciding node(s)
                for kc = 1:2
                    if onNode(kc)
                        coordno = nodesUC(curIx(kc));
                        coordOnNode(coordno) = fnodes(locpos(kc)/2);
                    end
                end
            elseif G.griddim == 2
                if round(tv) == 0
                    coordno = nodesUC(curIx);
                    coordOnNode(coordno) = fnodes(1);
                end
            end
        end
    elseif nn(k) >= 2
        unexpectedCount = unexpectedCount +1;
    end
end
if unexpectedCount > 0
    dispif(mrstVerbose, 'Ignoring %d face-splits due to unknown cause ...\n', ...
            unexpectedCount);
end

% if there are any coords on grid-nodes, we need to remove the corresponding
% poly-nodes and re-index poly-nodes
nodeReIndex = (1:size(coordOnNode))';
if ~all(isnan(coordOnNode))
    coix   = ~isnan(coordOnNode);
    nodeNo = coordOnNode(coix);
    % we could remove these coords, but don't bother now (stowaways)
    % use convention that negative node-no points to grid node, positive to
    % poly-node
    nodeReIndex(coix) = -nodeNo;
end
    
% re-index for cut-cells
cutCells.nodes = nodeReIndex(cutCells.nodes);

% make cut-faces
ix = cumsum([1;nn]);
ix = ix(isSplitFace);
if G.griddim == 3
    cutFaces = struct('ix',      fix(isSplitFace), ...
        'nodes',   nodeReIndex([nodesUC(ix), nodesUC(ix+1)]), ...
        'nodePos', [localNodePos(ix), localNodePos(ix+1)]);
elseif G.griddim == 2
    cutFaces = struct('ix',      fix(isSplitFace), ...
        'nodes',   nodeReIndex(nodesUC(ix)), ...
        'nodePos', localNodePos(ix));
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

%--------------------------------------------------------------------------
%{
function f = facesFromNeighbors(G, c)
if all(size(c)==[2,1])
    c = c';
end
nf = size(c,1);
f  = nan(nf,1);
fp = G.cells.facePos;
nFailed = 0;
for k =1:nf
    [c1, c2] = deal(c(k,1), c(k,2));
    if c1 == c2
        error('something wrong')
    else
        ix1 = fp(c1):(fp(c1+1)-1);
        ix2 = fp(c2):(fp(c2+1)-1);
        r = intersect(G.cells.faces(ix1), G.cells.faces(ix2));
        if ~(numel(r)==1)
            % disregard this situation, might be just a degenerate polygon
            % since validation is not done yet
            nFailed = nFailed +1;
        else
            f(k) = r;
        end
    end
end
if nFailed > 0
    dispif(mrstVerbose, 'Sub-function ''facesFromNeighbors'' in %s%s', ...
          '''extractGridInfo'' failed \n to identify %d faces.,', ...
          '(Not critical information)\n' , nFailed);
end
f = f(~isnan(f));
end
%}
%--------------------------------------------------------------------------
%{
Not currently used, but maybe useful

function g = createSliceGrid(G, cutCells, coords)
cix   = cutCells.ix;
nc    = numel(cix);
nn    = cutCells.nNodes;
nodes = cutCells.nodes;
% need to check for neagative node-indices for which we must collect coords
% from G
nix = nodes < 0;
if any(nix)
    [gnodes, ia, ib] = unique(-nodes(nix));
    ngNo = max(0, max(nodes)) + (1:numel(ia))';
    nodes(nix) = ngNo(ib);
    coords = [coords; G.nodes.coords(gnodes,:)];
end
% need segments
s1 = nodes;
pos = cumsum([1;nn(1:end)]);
s2  = getNext(s1, pos);
%s2 = [s1(2:end);nan];
%first = cumsum([1;nn(1:end-1)]);
%last  = cumsum(nn);
%s2(last) = s1(first);
segs = [s1, s2];
cellNo = rldecode((1:nc)', nn);
% cut from makePlanarGrid --------------------------------------------- 
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
%----------------------------------------------------------------------
g.nodes = struct('num', size(coords,1), 'coords', coords);
nf = size(edges,1);
g.faces = struct('num', nf, 'nodePos', (1:2:(2*nf+1))', 'neighbors', neigh, ...
                 'tag', zeros(nf, 1), 'nodes', reshape(edges', [], 1));
g.cells = struct('num', nc, 'facePos', cumsum([1; nn]), 'indexMap', nan(nc,1), ...
                 'parent', cix, 'faces', [edgenum, zeros(numel(edgenum),1)]);
g.griddim = 2;
g.type = {'sliceGrid'};
vst = mrstVerbose;
mrstVerbose(0);
g = computeGeometry(g);
g = repairNormals(g);
mrstVerbose(vst);
% normals need to be fixed!
g.type = [g.type, {'faulty face-normals'}]; 
g.faces.normals = nan;
% add cell normals
g.cells.normals = cutCells.normals;
end    
%}
%----------------------------------------------------------------------


