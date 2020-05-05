function fix = getFacesCloseToSegment(G, seg, varargin)
% Get index to faces that lie "in the vicinity" of the line segment seg. 
% 
% SYNOPSIS:
%  fix =  getFacesCloseToSegment(G, seg, 'pn1', pv1)
%
% DESCRIPTION:
% With default options, returned face-list will contain all but most likely 
% *not* only the faces actually intersected by the segment. Function is useful 
% for quickly obtaining a small subset for subsequent more accurate and costly
% intersecion evaluatioins.
% 
%
% REQUIRED PARAMETERS:
%
%   G      - Grid structure with geometry and preferably 'bbox'-field included
%            for G.faces (see addBoundingBoxFields)
%
%   seg    - Segment given as two points (2x3). If more than two points (nx3), 
%            line of best fit is computed with search radius expanded
%            according to maximal point-line distance. If point dimension
%            is 2, computations are perfomed in xy.
%
% OPTIONAL PARAMETERS:
%
% 'faceIx'  : consider only subset of grids faces (default all)
% 'fac'     : "Closeness"-factor. Default is 2/3 which is sufficient to 
%             gurantee that output includes all faces intersected by segment  
%
% RETURNS
%  fix      : resulting index of grid faces 
%
% SEE ALSO
%  computeTraversedCells, computeVerticalGridIntersection


opt = struct('faceIx',      [], ...
             'fac',        2/3, ...
             'tol',  sqrt(eps), ...
             'subsets',     []);

opt = merge_options(opt, varargin{:});
tol = opt.tol;

[dims, faceIx] = deal(':');

if size(seg,2) < 3
    dims = 1:2;
end
if ~isempty(opt.faceIx)
    faceIx = opt.faceIx;
end

if ~isfield(G.faces, 'bbox')
    warning('Precompute bounding box for faces for improved performance');
    G = addBoundingBoxFields(G);
end

bbox = G.faces.bbox(faceIx, dims);
cent = G.faces.centroids(faceIx, dims);


% if more than one segment, find a line of best fit
r = zeros(1,size(seg,2));
if size(seg,1) > 2
    p0 = mean(seg);
    pr = bsxfun(@minus, seg, p0);  %  relative coordinates
    [~, ~, v] = svd(pr, 0);        
    
    v = v(:,1)';
    % find required end-points and extra radius
    t    = pr*v';
    seg = [p0+min(t)*v; p0+max(t)*v];
    % distance points-to-line
    d = pr - t*v;
    r = max(abs(d)); % maximal distances
end

% main
% include points closer than 2/3 of bbox (1/2 is sufficient for cart grid)
fac  = opt.fac;
p1   = seg(1,:);
v    = diff(seg);

% all centroids relative to p1
m  = bsxfun(@minus, cent, p1);
% t-values for projections to line t*v 
t  = (m*v')/(v*v');
% extra length due to 2/3*bbox

% expand segment to include fac*bbox
dt = bbox*abs(fac*v')/norm(v);
% index to all centroids that project to segment
ix = find(t > -dt/norm(v)-tol & t < 1+dt/norm(v)+tol);

if isempty(ix)
    fix = [];
    return;
end

% do rest on subset ix
dv  = m(ix,:) - t(ix)*v;
% distance to segment
d   = sqrt(dot(dv, dv, 2));
% include bbox in dv-direction and point max distance to line
bbp = dot(bsxfun(@plus, fac*bbox(ix,:), r), abs(dv), 2);
ix2 = d.^2 <= bbp+tol; % squared since we want unit vector dv

% if we are (un)lucky, d is almost zero (line goes through centroid), include
% this situation

% dzero = d<sqrt(eps)*norm(seg);
% if any(dzero)
%     ix2 = ix2 | dzero;
% end

ix  = ix(ix2);
if ischar(faceIx)
    fix = ix;
else
    fix = faceIx(ix);
end
end
    