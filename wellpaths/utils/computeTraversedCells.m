function [T, extra] = computeTraversedCells(G, traj, varargin)
% Compute intersection details between grid and piecwise linear trajectory 
% 
% SYNOPSIS:
%   S = computeTraversedCells(G, traj, 'pn1', pv1)
%
% DESCRIPTION:
%
%
% REQUIRED PARAMETERS:
%
%   G      - Grid structure with geometry and preferably 'bbox'-field included
%            for G.faces (see addBoundingBoxFields)
%
%   traj   - piecewise trajectory given as list of points (nx3)
%
% OPTIONAL PARAMETERS:
%   faces  - Indices to subset of G.faces. Default is (1:G.faces.num)'
%
%   tol    - Relative tolerance for intersection testing. Should be > 0 to 
%            avoid effect of round-off.
%
%   exteriorFaceCorrection - if true, output T.weight is modified to account for 
%            exterior faces. This is done by computing the fracions of the
%            area-components normal to T.vec open to flow.
% 
% RETURNS:
%   T      - struct containing details of all traversed grid cells:
%            'cell'   - cell indices
%            'vec'    - vector decribing direction and length of cell-segment 
%            'weight' - weight of cell-segment, only ~=1 if segment is
%                       shared by more than one cell 
%
%  extra   - unprocessed output for debugging purposes etc 

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

opt = struct('faces',                           [], ...
             'tol',                      sqrt(eps), ...
             'exteriorFaceCorrection',       false, ...
             'segmentTol',               sqrt(eps), ...
             'removeDuplicates',             true, ...
             'duplicateTol',                  nan);
opt = merge_options(opt, varargin{:});

% remove part of traj outside grid bounding box
if ~isfield(G, 'bbox')
    G.bbox = [min(G.nodes.coords); max(G.nodes.coords)];
end
traj = trimTrajectory(traj, G.bbox);
        
nseg = size(traj,1)-1;

% get inital relevant subset of faces for whole trajectory
fix0 = getFacesCloseToSegment(G, traj, 'faceIx', opt.faces);
if isempty(fix0)
    error('Trajectory does not appear to intersect grid');
end

cnt       = 0;      
ec = {cell(100,1)};
T = struct('cell', ec, 'length', ec, 'p1', ec, 'p2', ec, 'weight', ec);

% Current segment
seg = struct('coord', [], 'curCell', 0, 'prevCell', 0, 'onFace', false);
for k = 1:nseg
    seg.coord = traj([k, k+1],:);
    while ~isempty(seg.coord)
        if any(seg.curCell > 0) && seg.onFace % entering new cell(s)
            cnt = cnt + 1;
            nc = numel(seg.curCell);
            T.cell{cnt}    = seg.curCell;
            T.p1{cnt}      = ones(nc,1)*seg.coord(1,:);
            T.p2{cnt}      = ones(nc,1)*seg.coord(2,:);
            T.weight{cnt}  = ones(nc,1)/nc;
        end
        if seg.curCell == 0     
            % not currently inside grid, get new face subset
            fix = getFacesCloseToSegment(G, seg.coord, 'faceIx', fix0);
            if isempty(fix)
                break;
                % might happen if segment is inside single cell
            end
        else
            % get all faces for current cell(s)
            fpos = G.cells.facePos(seg.curCell(1) + (0:1)');
            fix  = G.cells.faces(fpos(1):(fpos(2)-1));
            for cno = 2:numel(seg.curCell)
                fpos = G.cells.facePos(seg.curCell(cno) + (0:1)');
                fix = [fix, G.cells.faces(fpos(1):(fpos(2)-1))]; %#ok not likely
            end
        end
       
        % get face triangulation
        [v1, v2, v3, faceNo] = getTriangles(G, fix);
        
        % edges
        e1 = v2-v1;
        e2 = v3-v1;
        
        % some geometry (see moller-trumbore)
        t = seg.coord(1,:)-v1;
        d = repmat(diff(seg.coord), size(t, 1), 1);
        
        sg = cross(d, e2, 2);
        q = cross(t, e1, 2);
        den = dot(sg ,e1, 2);
        
        s = dot(q, e2, 2)./den;
        u = dot(sg, t,  2)./den;
        v = dot(q, d,  2)./den;
        
        % tolerance should be >0 to avoid missing intersections
        tol = opt.tol;
        % need 0 <= s <= 1 to be along segment
        ix = s >= -opt.segmentTol & s <= 1+opt.segmentTol;
        
        % intersection-test in barycentric coord
        ix = ix & u >= -tol & v >= -tol & (u+v) <= 1+tol;
        
        foundIntersection = false;
        
        if any(ix) % have hit one or more face
            ix = find(ix);
            [s, order] = sort(s(ix));
            ix = ix(order);
            
            % intersected faces
            f = faceNo(ix);
            
            % get neighbors in order according to direction of segment
            N = G.faces.neighbors(f,:);
            % Need to condsider normals of triangles. Assuming inner
            % product of these with face normal > 0 seems fair to get correct direction.
            nt = cross(e1(ix,:),e2(ix,:), 2);
            n  = G.faces.normals(f,:);
            nt = bsxfun(@times, sign(dot(n,nt,2)), nt);
            sgn = sign(nt*(diff(seg.coord)'));
            N(sgn==-1,:) = N(sgn==-1,[2 1]);
            
            % check if any intersection is outward
            if any(any(N(:,1)==seg.curCell(:).')) || all(seg.curCell == 0)
                foundIntersection = true;
            end
        end
        if ~foundIntersection % no new cell, add length to previous
            if seg.curCell > 0
                nc        = numel(seg.curCell);
                T.p2{cnt} = ones(nc,1)*seg.coord(2,:);
            end
            % we are done with segment
            seg.coord  = [];
            seg.onFace = false;
        else % deal with intersections
            if any(seg.curCell > 0)
                seg_new = findNextIntersection(seg, s, f, N, tol);
                % p_new still in same cell as p or on face of cell, add length
                if any(seg.curCell > 0)
                    nc         = numel(seg.curCell);
                    T.p2{cnt}  = ones(nc,1)*seg_new.coord(1,:);
                end
                seg = seg_new;
            else % deal with multiple intersections (we are currently 
                 % outside grid, or/and at k=1
                [seg, cells, x1, x2, w] = findMultipleIntersections(seg, s, f, N, tol);
                % lump together everything except last piece
                nEnd = numel(seg.curCell);
                if numel(cells) > 1
                    cnt = cnt + 1;
                    T.cell{cnt}    = cells(1:end-nEnd);
                    T.p1{cnt}      = x1(1:end-nEnd,:);
                    T.p2{cnt}      = x2(1:end-nEnd,:);
                    T.weight{cnt}  = w(1:end-nEnd,:);
                end
                if numel(cells) > 0
                    cnt = cnt +1;
                    T.cell{cnt}   = cells(end-nEnd+1:end);
                    T.p1{cnt}     = x1(end-nEnd+1:end, :);
                    T.p2{cnt}     = x2(end-nEnd+1:end, :);
                    T.weight{cnt} = w(end-nEnd+1:end,:);
                end
            end
        end
    end
end

% concatenate
flds = fieldnames(T);
for k = 1:numel(flds)
    T.(flds{k}) = vertcat(T.(flds{k}){:});
end
if nargout > 1
    extra = T;
end

% process/clean up cell segments
T = processCellSegments(T, tol);

if opt.exteriorFaceCorrection
    w = computeWICorrectionFactors(G, T.cell, T.vec);
    T.weight = T.weight.*w;
end
end

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

function [v1, v2, v3, faceNo] = getTriangles(G, f)
isSubset = true;
if isempty(f)
    f = (1:G.faces.num);
    isSubset = false;
end

if ~isSubset % full grid
    faceNo  = rldecode(1:G.faces.num, diff(G.faces.nodePos), 2) .';
    p       = G.faces.nodePos;
    next    = (2:size(G.faces.nodes, 1)+1) .';
    next(p(2 : end) - 1) = p(1 : end-1);
    
    v1 = G.nodes.coords(G.faces.nodes,:);
    v2 = G.nodes.coords(G.faces.nodes(next),:);
    v3 = G.faces.centroids(faceNo,:);
else
    [np1, np2] = deal(G.faces.nodePos(f), G.faces.nodePos(f+1));
    faceNo = rldecode(f(:), np2-np1);
    
    nodes = G.faces.nodes(mcolon(np1, np2-1));
    nnode = np2-np1;
    locpos = [1; cumsum(nnode)+1];
    next   = (2:numel(nodes)+1) .';
    next(locpos(2 : end) - 1) = locpos(1 : end-1);
    
    % triangle vertices
    v1 = G.nodes.coords(nodes,:);
    v2 = G.nodes.coords(nodes(next),:);
    v3 = G.faces.centroids(faceNo,:);
end
end

%--------------------------------------------------------------------------

function seg_new = findNextIntersection(seg, s, f, N, tol)
seg_new = seg; 
% intersection points pointing outwards
outIx  = any(N(:,1) == seg.curCell.', 2);
% select first occuring outIx in case of multiple intersections due to
% curved face
if nnz(outIx) > 1
    i1 = find(outIx, 1, 'first');
    ii = s-s(i1) <= tol;
    [s, f, N, outIx] = deal(s(ii), f(ii), N(ii,:), outIx(ii));
end
    
if ~any(outIx) % end of trajectory inside current cell
    seg_new.onFace = false;
    seg_new.coord = [];
else
    if ~all(outIx)
        [s, f, N] = deal(s(outIx), f(outIx), N(outIx,:));
    end
    ns = numel(s);
    ix = false(ns,1);
    smax = max(s);
    if numel(s) == 1 % single intersection, all well
        ix(1) = true;
    else % multiple intersections
        % get unique faces
        [f, fix] = unique(f, 'stable');
        if numel(f) < ns
            [s, N] = deal(s(fix), N(fix,:));
        end
        % select all points close to max(s) where seg is pointing outwards
        ix = abs(smax-s) < tol;
    end
   seg_new.curCell   = N(ix,2);
   seg_new.onFace = true;
   seg_new.coord(1,:) = seg.coord(1,:) +  smax*(diff(seg.coord));
end
end

%--------------------------------------------------------------------------

function [seg_new, cells, x1, x2, weight] = findMultipleIntersections(seg, s, f, N, tol)
ns = numel(s);
% unique faces (same face-index and almost same s)
six = cumsum([0; abs(s(2:end)-s(1:end-1))>tol/10]);
[f, fix] = unique([f, six], 'rows', 'stable');
f = f(:,1);
if numel(f) < ns
    [s, N] = deal(s(fix), N(fix,:));
    ns = numel(s);
end

if all(N(1:end-1, 2) == N(2:end, 1))
    % everything is in order :)
    s1   = [0; s];
    s2   = [s; 1];
    cells = [N(:,1); N(end,2)];
    weight = ones(ns+1, 1);
else % we have more troublesome intersections
    r = sortrows([ [N(:,1); N(:,2)],        ...
                   repmat((1:ns)', [2,1]),  ...
                   [true(ns,1); false(ns,1)] ]);
    [smin, smax] = deal(min(s), max(s));
    % loop through each cell present
    cnt = 0;
    w = zeros(ns+1, 1);
    s_ext = [0; s; 1];
    [cells, p1, p2] = deal(nan(size(r,1), 1));
    k=0;
    while k < size(r,1)
        k = k+1;
        cnt = cnt +1;
        k1  = k;
        cells(cnt) = r(k,1);
        while k < size(r,1) && r(k+1,1) == cells(cnt)
            k = k+1;
        end
        ix = r(k1:k, 2);
        sk = s(ix);
        if numel(ix)==1 % at start or end
            if (abs(sk-smin) <= tol) && r(k,3) % start
                pp = [0, 1] + 1;
            elseif (abs(sk-smax) <= tol) && ~r(k,3) % end
                pp = [ns, ns+1] +1;
            else
                error('Debug here')
            end
        else % multiple indices (two or more, select endpoints)
            [~, ixkMin] = min(sk);
            [~, ixkMax] = max(sk);
            pp = ix([ixkMin, ixkMax]) + 1;
        end
           % add weight to segments
           [p1(cnt), p2(cnt)] = deal(pp(1), pp(2));
           w(pp(1):(pp(2)-1)) =  w(pp(1):(pp(2)-1)) + 1;
    end
    [p1, p2, cells] = deal(p1(1:cnt), p2(1:cnt), cells(1:cnt));
    [s1, s2] = deal(s_ext(p1), s_ext(p2));
    weight = ones(cnt, 1);
    ds  = diff(s_ext);
    dsw = w.*ds;
    for k = 1:numel(cells)
        weight(k) = sum(dsw(p1(k):(p2(k)-1)))./sum(ds(p1(k):(p2(k)-1)));
    end
    
    weight = 1./weight;
    weight(~isfinite(weight)) = 1;
    % reorder according to segment starting points
    [s1, order] = sort(s1);
    [s2, cells] = deal(s2(order), cells(order));
end
% compute points
x1 = seg.coord(1,:) + s1*(diff(seg.coord));
x2 = seg.coord(1,:) + s2*(diff(seg.coord));
% update seg
seg_new = seg;
seg_new.curCell = cells(abs(s2-1)<tol);
seg_new.coord = [];
seg_new.onFace = false;
end
    
%--------------------------------------------------------------------------
function [T_out, pnt] = processCellSegments(T, tol)
% segment lengths
v    = T.p2-T.p1;
l    = sqrt(sum(v.^2,2));
ltol = tol*sum(l);
% remove short segments and outside grid
keep = T.cell > 0 & l >= ltol;
[c, v, w] = deal(T.cell(keep), v(keep,:), T.weight(keep));
% chech for repeated cells
[c, uix, ib] = unique(c, 'stable');
if numel(c) < numel(ib)
    % add vectors
    v = cellfun(@(x)accumarray(ib, x), mat2cell(v, size(v,1), [1,1,1]), ...
                'UniformOutput', false);
    v = horzcat(v{:});
    % average weights
    l  = l(keep);
    w  = accumarray(ib, w.*l)./accumarray(ib,l);
end
T_out = struct('cell', c, 'vec', v, 'weight', w);
if nargout > 1
    pnt = T.p1(uix,:);
end

end

%--------------------------------------------------------------------------

function w = computeWICorrectionFactors(G, cell, vec)
% Simple adjustmentfactor for external cell-faces.
% compute area components normal to vec and take fraction of interior to total 
w = ones(numel(cell), 1);
for k = 1:numel(cell)
    fpos = G.cells.facePos(cell(k) + (0:1)');
    fix  = G.cells.faces(fpos(1):(fpos(2)-1));
    N    = G.faces.neighbors(fix,:);
    if any(N(:)==0)
        %cf = G.faces.centroids(fix, :);
        vn = vec(k,:).'/norm(vec(k,:));
        % get area-weighted face normals        
        nf = G.faces.normals(fix, :);
        % projected normals
        nf = nf - (nf*vn)*vn';
        A  = sqrt(sum(nf.*nf,2));
        % interior faces     
        ie = all(N>0, 2);
        w(k) = sum(ie.*A)./sum(A);
    end
end
end

%--------------------------------------------------------------------------

function traj = trimTrajectory(traj, bbox)
% trim trajectory to to relevant region defined by bbox
% Include each trajectory point that is inside box or is part of segment that 
% intersects one of the box planes (this does not necessarily remove all 
% unwanted points, but this is just a rough initial trimming to get rid of 
% e.g., points along trajectory above formation)  
np = size(traj,1);
assert(np >= 2, 'Trajectory needs at least two points');
isBelow = bsxfun(@le, traj, bbox(1,:));
isAbove = bsxfun(@ge, traj, bbox(2,:));
isInside   = all( ~isBelow & ~isAbove ,2);
isCrossing = any(xor(isAbove(1:np-1,:),isAbove(2:np,:)) | ...
                 xor(isBelow(1:np-1,:),isBelow(2:np,:)), 2);
include = isInside | [false; isCrossing] | [isCrossing; false];

if ~any(include)
    error('Trajectory does not intersect the bounding box of the grid');
end
traj = traj(include, :);
end
