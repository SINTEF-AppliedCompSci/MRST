function [vertices,fractures] = mergeFractures(vertices,fractures,vert_new,frac_new,box,precision)
% merge two set of fractures. To avoid silver cells
% new points close to old fractures are projected onto the old fractures
% if both endpoint of a fracture is snapped to the same old fracture
% it is removed to avoid duplicated fractures.
%
% SYNAPSES
% [vertices,fractures] = mergeFractures(vertices,fractures,vert_new,frac_new,box,precision)
%
% PARAMETERS:
%
%       vertices  - original set of vertices
%       fractures - original set of fractures
%       vert_new  - new set of vertices
%       frac_new  - new set of fractures
%       box       - bounding box of the domian
%       presicion - the releative presicion of the grid
%
% OUTPUT:
%       vertices  - the merged vertices
%       fractures - the merged fractures
%
% Copyright 2011-2012 University of Bergen, 2013 IRIS AS
%
% This file is licensed under the GNU General Public License v3.0.


% do nothing if their is no old fractures
if isempty(vertices)
    vertices = vert_new;
    fractures = frac_new;
    return
end

% snap close child vertices to parent edges.
[dist,ind] = distP2E(vert_new,vertices,fractures);

% dist_tol is defined relative to domain size
dist_tol = 2 * norm(diff(box)*precision);

% prodject all nodes closer than the given tolerance.
isprodnodes = dist<dist_tol;
prodnodes =find(isprodnodes);

% do the prodjection
pv = zeros(numel(prodnodes),2);
% prodjekt point to its closest edge
for i = 1 : numel(prodnodes)
    edge = fractures(ind(prodnodes(i)),:);
    pv_this = projectP2Line(vertices(edge(1),:),vertices(edge(2),:),vert_new(prodnodes(i),:));
    pv_this = snap_to_line(pv_this,vertices(edge(1),:),vertices(edge(2),:),box,struct('precision', precision));
    pv(i,:) = pv_this;

end

% remove all fractures where both points are snapped to
% the same edges as they are now parallel to their edges
% in this way we may lose track of tags
% TODO: duplicate tags to make sure that all tags are kept
remed = ind(frac_new(:,1)) == ind(frac_new(:,2)) & all(isprodnodes(frac_new(:,1:2)),2);
frac_new = frac_new(~remed,:);
vert_new(isprodnodes,:) = pv;

% add the new fractures and vertices to the old ones.
numVertices = size(vertices,1);
vertices = [vertices; vert_new];
frac_new(:,1:2) = frac_new(:,1:2) + numVertices;
fractures = [fractures; frac_new];

end

%% helpers
function  [dist,ind] = distP2E(points,vertices,edges)
% A wrapper around distance_to_closest_line.m

x0 = vertices(edges(:,1),1);
y0 = vertices(edges(:,1),2);
x1 = vertices(edges(:,2),1);
y1 = vertices(edges(:,2),2);

% call mex function
[dist,ind] = distance_to_closest_line(points(:,1),points(:,2),x0,y0,x1,y1);

% alternative .m function to use (vectorized)
% d_new = distance_to_line_segments2(old_points(:,1),old_points(:,2),x0,y0,x1,y1);
% [dist,ind] = min(d_new,[],2);

% % find the nearest distance from the points to the lines given by
%
% dist = inf(size(points,1),1);
% ind = zeros(size(points,1),1);
% % Loop over all lines
%   for i = 1:size(edges,1)
%       % Find the distance from all points to the line
%
%       d_new = distance_from_line_segment(points,vertices(edges(i, 1:2), :));
%
%       % Have we found a new minimum distance?
%       [dist,I] = min([dist,d_new],[],2);
%       ind(I==2) = i;
%
%   end
end

function p = projectP2Line(here,there,p)
% project the point to its edge
v = there-here;
c = p-here;
p = bsxfun(@plus,bsxfun(@times,dot(c,v,2)./dot(v,v,2),v),here);
end

function [vertices] = snap_to_line (vertices, p1,p2, box,opts)
% Move vertices to the closest (structured) grid point. Use this function to
% avoid having more than one point within each (fine) grid cell.  %edge = sort(edge,2);
precision = sqrt(2) * norm(diff(box))*opts.precision/norm(p2-p1);

% dist_tol is defined relative to domain size
if precision > eps && ~isempty (vertices),
    a = round(norm(vertices-p1)/norm(p2-p1)/precision)*precision;
    a = max(a,precision);
    a = min(a,1-precision);
    vertices = a * (p2-p1)+p1;

end

end
