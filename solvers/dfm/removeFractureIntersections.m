function [vertices, edges] = removeFractureIntersections(vertices, edges, box, opts)
% Partition intersecting edges by adding a point in the intersection.
%
% Copyright 2011-2012 University of Bergen
% This file is licensed under the GNU General Public License v3.0.

% Column index of x and y coordinate in vertices
x = 1; y = 2;
% Loop over all edges, search for intersections. The number of edges can
% change, thus we use while
i = 1;

while i <= size(edges,1)

    % To test if two line segments intersect is somewhat expensive.
    % Therefore, we take a two-stage approach: First, do a coarse search to
    % disregard edges that clearly do not intersect edge i. For all
    % remaining edges (usually not that many), we test more carefully if
    % there is a crossing.
    % The coarse search is divided into two parts: Firstly, edges that do not
    % cross the extension of edge i are removed. Secondly, we remove edges
    % where the circurmscribed circles are not overlapping.
    % These two algorihms remove most of the edges, and they can easily be
    % vectorized. Other criteria for sorting may be more appropriate, who
    % knows.

    %% First search for line intersections

    % Find endpoints of all edges (subscript 1 and 2 refers to column
    % position in edges)
    px_1 = vertices(edges(:,1),x);
    py_1 = vertices(edges(:,1),y);
    px_2 = vertices(edges(:,2),x);
    py_2 = vertices(edges(:,2),y);

    % Coefficients needed to represent the line on the form ax + by = c
    a =   py_2 - py_1;
    b = -(px_2 - px_1);

    % Midpoints of line i
    xm = 0.5 * (px_1(i) + px_2(i));
    ym = 0.5 * (py_1(i) + py_2(i));

    %
    % For all lines, find which side of line i it's two endpoints are. If
    % c1 and c2 have different signs, they will be on different sides of
    % line i. See
    %
    % http://stackoverflow.com/questions/385305/efficient-maths-algorithm-to-calculate-intersections
    %
    % for more information. This works as a coarse search for intersection
    % points, since we consider i a line and not a line segment, se below.
    c1 = a(i)*(px_1 - xm) + b(i)*(py_1 - ym);
    c2 = a(i)*(px_2 - xm) + b(i)*(py_2 - ym);

    % Find all lines that intersect with the extension of line i.
    % Intersections will contain
    % both lines that are crossing (sign(c1)*sign(c2) < 0), and lines which
    % starts or end on line i. Intuitively, it seems that by including the
    % latter case, the robustness of the algorithm is increased (I'm not
    % sure of how susceptible this gridding algorithm is to rounding
    % errors), but it might lead to unnecessary operations as wel.
    lineIntersections = (sign(c1) ~= sign(c2));

    intersections = find(lineIntersections);

    counter = 1;
    while numel(intersections) > 0 && counter <= numel(intersections)

        % Process the first line that intersect with i
        j = intersections(counter);

        % Line j is an intersection if it crosses the extension of line i
        % (ie it crosses the infinite line that goes through the endpoints
        % of line i), but we don't know if it actually crosses the line
        % segment i. Now we do a more refined search to find if the line
        % segments intersects. Note that there is no help in vectorizing
        % lines_intersect and computing intersection points for all lines
        % in intersections, since line i may be split, and the
        % intersection points must recalculated. It may be possible to
        % reorganize this while-loop by computing all intersection points
        % (vectorized), and only recompuing if line i is split, but we
        % keep things simple for now.
        pt = lines_intersect(px_1(i),py_1(i),px_2(i),py_2(i),px_1(j),py_1(j),px_2(j),py_2(j));

        % If we have found an intersection point, the lines must be split,
        % and a new vertex might be introduced. There are two possible
        % cases here: One of the lines starts on the other line (eg i
        % represents a boundary and j is an interior constraint tha goes
        % all the way to the boundary, this can happen in the interior as
        % well). In this case, one edge is split, and no new vertices
        % are needed. The other option is we've found a new crossing, and
        % we must introduce a new vertex and split both edges.
        if numel(pt) > 0

            % Possibly split edge i
            [vertices, edges, delta_j_for_i] = split_edge (vertices, edges, i, pt, box, opts);

            % If edge i was split, a new edge was introduced in line i + 1.
            % Since j > i, we may have to increase j here.
            j = j + delta_j_for_i;

            % Possibly split edge j
            [vertices, edges, delta_j_for_j] = split_edge (vertices, edges, j , pt, box, opts);

            % The introduction of new edges have modified edge indices, thus
            % intersections must be updated.
            intersections = intersections + delta_j_for_i + delta_j_for_j;

            % If any of the lines were split, new points might be
            % introduced, and we need to update pxy_12
            if delta_j_for_i + delta_j_for_j > 0
                px_1 = vertices(edges(:,1),x);
                py_1 = vertices(edges(:,1),y);
                px_2 = vertices(edges(:,2),x);
                py_2 = vertices(edges(:,2),y);
            end

        end
        % Continue with the next element in intersections
        counter = counter + 1;
    end

    % We are done with all intersections of line i, move on to the next one
    i = i + 1;
end
