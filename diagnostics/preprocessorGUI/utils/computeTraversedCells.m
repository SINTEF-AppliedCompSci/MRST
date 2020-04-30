function tt = computeTraversedCells(G, traj, varargin)
opt = struct('faces',               [], ...
    'tol',                  0, ...
    'segmentTol',           0, ...
    'removeDuplicates',  true, ...
    'duplicateTol',       nan);
opt = merge_options(opt, varargin{:});

if isnan(opt.duplicateTol)
    opt.duplicateTol = opt.tol;
end

nseg = size(traj,1)-1;

[fix0, subsets] = getFacesCloseToSegment(G, traj, 'faceIx', opt.faces);

enterCell    = 0;
nextCellMode = false;
count        = 0;

ec = {cell(100,1)};
tt = struct('cells', ec, 'lengths', ec, 'coord1', ec, 'coord2', ec, 'weights', ec);
for k = 1:nseg
    curseg = traj([k, k+1],:);
    while ~isnan(curseg(1,1))
        count = count + 1;
        if all(enterCell == 0) || ~nextCellMode
            fix = getFacesCloseToSegment(G, curseg, 'faceIx', fix0, 'subsets', subsets);
            nextCellMode = false;
        else
            fpos = G.cells.facePos([enterCell;enterCell+1]);
            fix  = G.cells.faces(fpos(1):fpos(2));
            nextCellMode = true;
        end
       
        % get face triangulation
        [v1, v2, v3, faceNo] = getTriangles(G, fix);
        
        % edges
        e1 = v2-v1;
        e2 = v3-v1;
        
        % some geometry (see moller-trumbore)
        t = curseg(1,:)-v1;
        d = repmat(diff(curseg), size(t, 1), 1);
        
        p = cross(d, e2);
        q = cross(t, e1);
        den = dot(p ,e1, 2);
        
        s = dot(q, e2, 2)./den;
        u = dot(p, t,  2)./den;
        v = dot(q, d,  2)./den;
        
        % may use some tollerance to avoid missing intersections between faces
        tol = opt.tol;
        % need 0 <= s <= 1 to be along segment
        ix = s >= -opt.segmentTol & s <= 1+opt.segmentTol;
        
        % intersection-test in barycentric coord
        ix = ix & u >= -tol & v >= -tol & (u+v) <= 1+tol;
        
        % empty output if no intersection is found
        %[c, len, p1, p2, f, weight, startIx, endIx] = deal([]);
        if any(ix)
            ix = find(ix);
            [s, order] = sort(s(ix));
            ix = ix(order);
            %[s, u, v] = deal(s(ix), u(ix), v(ix));
            
            % intersected faces
            f = faceNo(ix);
            %f = f(order);
            
            % get neighbors in order according to direction of segment
            N = G.faces.neighbors(f,:);
            sgn = sign(G.faces.normals(f,:)*(d(1,:)'));
            N(sgn==-1,:) = N(sgn==-1,[2 1]);
            
            if opt.removeDuplicates
                % remove duplicate intersections due to neighboring triangles
                [uix, foundDuplicates] = getUniqueFaceIx(s, f, N, opt.duplicateTol);
                if foundDuplicates
                    [s, f, N, ix] = deal(s(uix), f(uix), N(uix,:), ix(uix));
                end
            end
            
            % get barycentric and cart coord
            [u, v] = deal(u(ix), v(ix));
            coord  = (1-u-v).*v1(ix,:) + u.*v2(ix,:) + v.*v3(ix,:);
            
            info = processCellSegments(s, N);
            %
            c = info.cells;
            % add endpoints
            coord = [curseg(1,:); coord; curseg(2,:)]; %#ok
            [p1, p2] = deal(coord(info.ix1+1,:), coord(info.ix2+1,:));
            weight = info.weights;
            [startIx, endIx] = deal(info.startCell, info.endCell);
            
            % also compute cellsegment lengths
            s_ext = [0; s; 1];
            len = (s_ext(info.ix2+1) - s_ext(info.ix1+1))*norm(diff(curseg));
            %    [c, p1, p2, f, weight, startCell, endCell]
            %    nseg = numel(s)-1;
            %    c    = N(1:nseg,2);
            %    p1   = coord(1:nseg, :);
            %    p2   = coord(2:nseg+1, :);
            
            % add segment first and last if nonempty
            %    if s(1) > opt.segmentTol
            %        c  = [N(1,1); c];
            %        p1 = [seg(1,:);   p1];
            %        p2 = [coord(1,:); p2];
            %    end
            %    if s(end) <= 1 - opt.segmentTol
            %        c  = [c; N(end, 2)];
            %        p1 = [p1; coord(end,:)];
            %        p2 = [p2; seg(2,:)];
            
            keep = len > sqrt(eps);
            if nextCellMode
                if c(info.endCell) == enterCell % still inside same cell
                    curseg = nan;
                else % gone through cell
                    if ~isempty(info.endCell)
                        keep(info.endCell) = false;
                        curseg(1,:) = p1(info.endCell(1), :);
                        enterCell = c(info.endCell);
                    else % outide grid
                        curseg(1,:) = p1(end, :);
                        enterCell = 0;
                    end
                end
            else % expect multiple and that all neccessary faces are present 
                curseg = nan;
                if ~isempty(info.endCell)
                    enterCell = c(info.endCell);
                else
                    enterCell = 0;
                end
                nextCellMode = true;
            end
            [c, len, p1, p2, weight] = deal(c(keep), len(keep), p1(keep,:), p2(keep,:), weight(keep));
            if count > 1 && tt.cells{count-1}(end) == c(1) % add length to previous
                [l1, l2] = deal(tt.lengths{count-1}(end), len(1));
                [w1, w2] = deal(tt.weights{count-1}(end), weight(1));
                tt.lengths{count-1}(end)    = l1+l2;
                tt.weights{count-1}(end)   = (l1*w1 + l2*w2)/(l1 +l2);
                tt.coord2{count-1}(end,:)  = p2(1,:); 
                ix = 2:numel(c);
                [c, len, p1, p2, weight] = deal(c(ix), len(ix), p1(ix,:), p2(ix,:), weight(ix));
            end
            if ~isempty(c)
                tt.cells{count}   = c;
                tt.lengths{count} = len;
                tt.coord1{count}  = p1;
                tt.coord2{count}  = p2;
                tt.weights{count}  = weight;
            else
                count = count-1;
            end
        else
            if all(enterCell > 0) && count > 1
                tt.lengths{count-1}(end) = tt.lengths{count-1}(end) + norm(diff(curseg));
                tt.coord2{count-1}(end,:)  = curseg(2,:);
            end
            count = count-1;
            %if nextCellMode % something whent wrong, redo
            %    nextCellMode = false;
            %    disp('obs')
            %else
            curseg = nan;
            %end
        end
    end
end
flds = fieldnames(tt);
for k = 1:numel(flds)
    tt.(flds{k}) = vertcat(tt.(flds{k}){:});
end
end
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

function [uix, flag] = getUniqueFaceIx(s, f, N, tol)
flag = false;
uix  = true(numel(s), 1);
% find potential clusters
d   = s(2:end)-s(1:(end-1));
cix = cumsum([1; double(d>tol)]);
if cix(end) < numel(s)
    for k = 1:cix(end)
        ii = cix==k;
        if nnz(ii) > 1
            [~, ia, ib] = unique([f(ii), N(ii,:)], 'rows');
            if numel(ia) < numel(ib)
                tmp     = false(numel(ib), 1);
                tmp(ia) = true;
                uix(ii) = tmp;
                flag = true;
            end
        end
    end
end
end

function info = processCellSegments(s, N)
% remap to local
cells = reshape(unique(N), [], 1);
if ~cells(1) == 0; cells = [0; cells]; end
remap = sparse(cells+1, 1, (1:numel(cells))');
N = full(remap(N+1));
N = N-1;
if size(N,2) == 1
    N = N';
end
if cells(1) == 0, cells = cells(2:end); end
nc = numel(cells);
ns = numel(s);
% consider line segments [start, s(1)], [s(2), s(2)], ..., [s(m) end]
% enter and exit points for each cell
cellPnts = [zeros(nc,1), (ns+1)*ones(nc,1)];
[nz1, nz2] = deal(find(N(:,1)>0), find(N(:,2)>0));
% assume at most one entry and exit per cell, if multiple are occuring
% select first entrypoint and last exitpoint
cellPnts(N(nz1,1),2) = nz1;
cellPnts(flipud(N(nz2,2)),1) = flipud(nz2);
% identical points may mess up the order, make sure increasing index
switchIx = cellPnts(:,1)>cellPnts(:,2);
cellPnts(switchIx,:) = cellPnts(switchIx, [2 1]);
% segments included in each cell
%cellSegs = [cellPnts(:,1)+1, cellPnts(:,3)];
segNo  = mcolon(cellPnts(:,1)+1, cellPnts(:,2))';
cellNo = rldecode((1:nc)', diff(cellPnts,[],2));
% number of cells for each segment:
cellCnt   = accumarray(segNo, ones(size(segNo)), [ns+1, 1]);
s_ext     = [0; s(:); 1];
segLen    = diff(s_ext);
segWeight = segLen./cellCnt;
segWeight(cellCnt==0) = 0;

cellWeight = accumarray(cellNo, segWeight(segNo))./accumarray(cellNo, segLen(segNo));
cellWeight(~isfinite(cellWeight)) = 0;

% finally order cells according to s-value
cellcent = ( s_ext(cellPnts(:,1)+1) + s_ext(cellPnts(:,2)+1) )/2;
[~, order] = sort(cellcent);

info = struct('cells', cells(order), ...
    'ix1',   cellPnts(order,1), ...
    'ix2', cellPnts(order,2), ...
    'weights', cellWeight(order), ...
    'startCell', find(cellPnts(order,1)==0), ...
    'endCell', find(cellPnts(order,2)==ns+1) );
% if isempty(info.startCell)
%     info.startCell = 0;
% end
% if isempty(info.endCell)
%     info.endCell = 0;
% end
end









%function cix = cluster(s, tol)
%d   = s(2:end)-s(1:(end-1));
%cix = cumsum([1; double(d>tol)]);
%end

