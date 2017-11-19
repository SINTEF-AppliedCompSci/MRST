function [vertices, edges, new_line] = ...
    split_edge(vertices, edges, i, pt, box, opts)
% Split an edge into two parts
%
% Create two new edges from one old. Split an edge at index i into two
% at the point (x, y) given by pt and insert the two parts in its
% place. (The rest of the edge set is adjusted . It assumes that the point
% (x, y) is on the line. The result variable new_line is 0 if no line is
% added since the point is already one of the endpoints of the old line, or
% 1 if a new line is created, i.e. the old line fissioned into two lines
% (the value is the number of lines that are added).
%
% Copyright (C) 2006-2007 Uni Research AS
% This file is licensed under the GNU General Public License v3.0.

% get the indices of the two points in the old edge
a = edges(i, 1);
b = edges(i, 2);

% each of the two new edges will be of the same type as the old one
t = edges(i, 3:end);

% make sure that we have a point for the new middle point, c
[vertices, c] = add_point (vertices, pt, box, opts);

% only split the line if we actually got a new point, since otherwise we
% would get an empty line. note that this code depends on the property of
% the add_point function to identify the same point for (approximately)
% the same coordinates
if ~((a == c) || (b == c))

    % insert these two new edges in the place where the old one used to be
    edges = [edges(1:i-1, :); ...
        a c t; ...
        c b t; ...
        edges(i+1:size(edges, 1), :)];
    new_line = 1;
else
    new_line = 0;
end;
