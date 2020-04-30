function [fix, subsets] = getFacesCloseToSegment(G, seg, varargin)
opt = struct('faceIx',            [], ...
             'crossSection',    'xy', ...
             'fac',              2/3, ...
             'subsets',           []);
opt = merge_options(opt, varargin{:});

[gddim, segdim] = deal(G.griddim, size(seg,2));

if isempty(opt.faceIx)
    faceIx = (1:G.faces.num)';
else
    faceIx = opt.faceIx;
end
if ~isfield(G.faces, 'bbox')
    warning('Precompute bounding box for faces for improved performance');
    G = addBoundingBoxFields(G);
end

% if more than one segment, find a representative line
r = zeros(1, segdim);
if size(seg,1) > 2
    p0 = mean(seg);
    pr = bsxfun(@minus, seg, p0);
    [~, ~, v] = svd(pr, 0);
    
    v = v(:,1)';
    % find required end-points and extra radius
    t    = pr*v';
    % if t(1) > 0, vector v is reverse direction of p1, p2, ..., pn
    %if t(1) > 0
    %    [t, v] = deal(-t, -v);
    %end
    seg = [p0+min(t)*v; p0+max(t)*v];
    d = pr + t*v;
    r = max(abs(d));
end

% main
% include points closer than 2/3 of bbox (1/2 is sufficient for cart grid)
fac  = opt.fac;
p1   = seg(1,:);
v    = diff(seg);

if ~isempty(opt.subsets)
    bbox = opt.subsets.bbox;
    cent = opt.subsets.centroids;
else
    bbox = G.faces.bbox(faceIx,1:segdim);
    cent = G.faces.centroids(faceIx,1:segdim);
end


m  = bsxfun(@minus, p1, cent);
t  = -(m*v')/(v*v');
% extra length due to 2/3*bbox
vn = v/norm(v);
dt = bbox*abs(fac*vn')/norm(v);
ix = find(t>-dt & t < 1+dt);
% do rest on subset ix
dv  = m(ix,:) + t(ix)*v;
d   = sqrt(dot(dv, dv, 2));
bbp = dot(bsxfun(@plus, bbox(ix,:), r/fac), abs(dv), 2)./d;
ix2 = d<=fac*bbp;

ix  = ix(ix2);
fix = faceIx(ix);
if nargout == 2
    subsets = struct('bbox', bbox(ix,:), 'centroids', cent(ix,:));
end
end
    