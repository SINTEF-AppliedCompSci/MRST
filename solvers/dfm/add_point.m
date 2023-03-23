function [vertices, i, new_point] = add_point(vertices, pt, box, opts)
% Add a new point to the cloud.
%
% Ensure the existence of a point in the vertex set if it is not already
% there. When the function return, the point will be present at
% vertices(i). new_point contains 0 if an old point was found, and 1 if a
% new point is created, i.e. it is the number of points added.
%
% Portions Copyright (C) 2006-2007 Uni Research AS
% This file is licensed under the GNU General Public License v3.0.


% coordinates are specified in floating-point values which by nature are
% inexact. we don't mind if we hit the points exactly; we only have to
% get close.
vertices = snap_to_grid (vertices, box, opts);
pt = snap_to_grid (pt, box, opts);

% check against each existing point. if all components in the row is equal,
% then return this index (there should be only one)
x = 1; y = 2;
if ~isempty (vertices),
    i = find (and (vertices(:, x) == pt(x), vertices(:, y) == pt(y)));
    % in case more than one vertex was found (since we did some coalescing in
    % the calls to snap_to_grid above), then only return the first
    if length (i) > 1, i = i(1); end;
else
    i = [];
end;

% if we don't find the point anywhere, then we'll add it to the end
if isempty (i),
    i = size (vertices, 1) + 1;
    vertices(i, :) = pt;
    new_point = 1;
else
    new_point = 0;
end;