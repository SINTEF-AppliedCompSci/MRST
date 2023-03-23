function [new_points] = ...
    remove_closepoints(vertices, edges, old_points, tol)
% Remove all points that are too close to the lines, since these will only
% cause slivers.

% SYNOPSIS
% [new_points] = ...
%       remove_closepoints (vertices, edges, old_points, space, opts)
%
% PARAMETERS
%   vertices    - a list of coordinates that makes the edges. Only 2d.
%   edges       - a list of lines/edges
%
%   old_points  - list of points that you want to remove close points from
%   tol         - either a scalar giving the tolerance
%                 or a vector with a tolerance value for each of the tag
%                 given in edges(:,3).
%
%  RETURNS
%
%   new_points  - a new list of points where close points are removed
%
% Copyright 2011-2012 University of Bergen
%
% This file is licensed under the GNU General Public License v3.0.

% get end points of the edges
x0 = vertices(edges(:,1),1);
y0 = vertices(edges(:,1),2);
x1 = vertices(edges(:,2),1);
y1 = vertices(edges(:,2),2);

% call mex function to find closest distance
[dist,ind] = distance_to_closest_line(old_points(:,1),old_points(:,2),x0,y0,x1,y1);


numPoints = size(old_points,1);

% either use the same tolerance for all edges
if numel(tol) == 1 || size(edges,2) < 3
    TOL = tol * ones(numPoints,1);
else % or use different tolerances for different tags
    tags = edges(ind,3);
    TOL = zeros(numPoints,1);
    for i = 1:length(tol);
        TOL(tags==i) = tol(i);
    end
    % make sure we have not forgot anyone
    assert(all(TOL)>0);
end

% Find points that are far enough away from all lines
hit = dist > TOL;

% only keep the points that are not too close
new_points = old_points(hit,:);


