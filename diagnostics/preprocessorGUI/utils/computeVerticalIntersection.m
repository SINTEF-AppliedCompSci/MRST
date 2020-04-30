function polygons = computeVerticalIntersection(G, seg)

seg  = seg(:, 1:2);
nseg = size(seg,1)-1;
% pick intial candidate subset of faces
if nseg > 1
    [fix0, subsets] = getFacesCloseToSegment(G, seg);
else
    [fix0, subsets] = deal([]);
end

[p, t, fno] = deal(cell(1, nseg));
ts = 0;

cumlen = zeros(nseg+1, 1);
for k = 1:nseg
    curseg = seg([k, k+1],:);
    fix = getFacesCloseToSegment(G, curseg, 'faceIx', fix0, 'subsets', subsets);
    
    [p1, p2, faceNo] = getFaceSegments(G, fix);
    
    % current vertical plane through curseg n*x+d=0;
    v = diff(curseg);
    n = [-v(2), v(1), 0];
    d = curseg(1,:)*n(1:2)';
    [pcur, ix] = getPlaneIntersections(p1, p2, n, d, curseg);
    % also compute length along plane
    rp = bsxfun(@minus, pcur(:, 1:2), curseg(1,:));
    tcur = ts + sqrt(dot(rp, rp, 2));
    ts = ts + norm(v);
    cumlen(k+1) = ts;
    [p{k}, t{k}, fno{k}] = deal(pcur, tcur, faceNo(ix));
end
[p, t, fno] = deal(vertcat(p{:}), vertcat(t{:}), vertcat(fno{:}));
% extend to cells
cno  = G.faces.neighbors(fno,:);
[nz1, nz2] = deal(cno(:,1)>0, cno(:,2)>0);
[p, t, cno] = deal([p(nz1,:); p(nz2,:)], [t(nz1); t(nz2)], [cno(nz1,1); cno(nz2,2)]);

% reindex c 1:n
cix = unique(cno);
if cix(1) == 0
    cix = cix(2:end);
end
rix = sparse(G.cells.num, 1);
rix(cix) = (1:numel(cix));
cno_loc = full(rix(cno));
npnt = accumarray(cno_loc, ones(size(cno_loc)));

npol = numel(npnt);
% compute angle in (t, p(:,3)) -plane of points wrt face center
pp = [t, p(:,3)];
cent  = [accumarray(cno_loc, pp(:,1)), accumarray(cno_loc, pp(:,2))]./npnt;
rp    = pp - cent(cno_loc,:);
theta = atan2(rp(:,2), rp(:,1));
[~, order] = sortrows([cno_loc, theta]);
% each point now occurs twice, will come in pairs after sorting
% here we will loose some poly neighbour-info, redo if this is important
%order = order(1:2:end); 
%npnt = npnt/2;
[p, t, cno_loc] = deal(p(order,:), t(order), cno_loc(order));
cno = cix(cno_loc);

% make vertice info
maxp = max(npnt);
unique(npnt)

cumc  = cumsum(npnt);
pos   = [1; cumc(1:end-1)+1];
nodes = zeros(npol, maxp);
cl    = (1:numel(cno))';
for k = 1:maxp
   nodes(:,k) =  cl(min(pos+k-1, cumc));
end

polygons = struct('nodes',           nodes, ...
                  'coords3D',            p, ... 
                  'coords2D',  [t, p(:,3)], ...
                  'cellIx',            cix, ...
                  'segments',          seg, ...
                  'cumlength',      cumlen);
end
% -------------------------------------------------------------------------
function [p1, p2, faceNo] = getFaceSegments(G, f)
[np1, np2] = deal(G.faces.nodePos(f), G.faces.nodePos(f+1));
faceNo = rldecode(f(:), np2-np1);

nodes = G.faces.nodes(mcolon(np1, np2-1));
nnode = np2-np1;
locpos = [1; cumsum(nnode)+1];
next   = (2:numel(nodes)+1) .';
next(locpos(2 : end) - 1) = locpos(1 : end-1);

% triangle vertices
p1 = G.nodes.coords(nodes,:);
p2 = G.nodes.coords(nodes(next),:);
end

function [p, ix] = getPlaneIntersections(p1, p2, n, d, seg)
v  = p2-p1;
vn = v*n';
ok = abs(vn) >0;
t  = nan(size(vn));
%t  = (p1*n'-d)./v*n';
t(ok) = -(p1(ok,:)*n'-d)./vn(ok);
ix = find(t>=0 & t <= 1);
p  = p1(ix,:) + bsxfun(@times, t(ix), v(ix,:));
% only those along segment
v = diff(seg);
s = (p(:,1:2)*v(1:2)' - seg(1,:)*v')/(v*v');
ix2 = s>=0 & s <=1;
[p, ix]  =deal(p(ix2, :), ix(ix2));
end