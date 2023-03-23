function [vertices, new_edges, map] = ...
                        partition_edges(vertices, old_edges, space, box, opts)
% [vertices, new_edges, num_added] = ...
%                      partition_edges (vertices, old_edges, space, box, opts)
%
% Make sure that no edges are longer than the give size. If an edge is
% longer, it is split up into subsegments which are then joined together.
% map is a m-by-2 array which contains in the first column the index of the
% vertex added and in the second column the index of a line to which it is
% constrained.
%
% Portions Copyright (C) 2006-2007 Uni Research AS
% This file is licensed under the GNU General Public License v3.0.

   new_edges = [];
   map = [];

   % do for each edge in the edge set
   for i = 1:size (old_edges, 1)
      % get the coordinates of each endpoint of the line
      x_a = vertices(old_edges(i, 1), 1);
      y_a = vertices(old_edges(i, 1), 2);
      x_b = vertices(old_edges(i, 2), 1);
      y_b = vertices(old_edges(i, 2), 2);

      % each partition is going to be of the same type as the main edge
      c = old_edges(i, 3:end);

      % figure out the number of points we'll need in addition to the
      % starting point (the endpoint counts, so there should always be at
      % least one) to describe this line in terms of smaller segments with
      % at most 'space' length
      total_length = sqrt ((x_b - x_a)^2 + (y_b - y_a)^2);
      num_points = ceil (total_length / space);

      % parameterize the line so that it is described according to a
      % counter going from 0 to num_points:
      %     x = t/n * (x_b - x_a) + x_a
      %     y = t/n * (y_b - y_a) + y_a
      dx = (x_b - x_a) / num_points;
      dy = (y_b - y_a) / num_points;

      % start out at the first point of the old line
      last_point = old_edges(i, 1);

      % create new points and setup segments between them for all interior
      % points on the line
      for t = 1:num_points-1
         % calculate the coordinates of the point according to the
         % parameterized formula described above
         x = dx * t + x_a;
         y = dy * t + y_a;

         % make sure that we actually have one such point present; since
         % the edges are assumed to be without crossings at this point,
         % there shouldn't be any points on the line we could possible us;
         % this call should always generate a new point
         [vertices, new_point] = add_point (vertices, [x y], box, opts);

         % draw another segment of the line, and let the pen be ready at
         % the ending point of this segment for the next.
         new_edges = [new_edges; last_point, new_point, c];
         last_point = new_point;

         % the line to which the point is constrained will always be the
         % last line we have added so far, and the index is in the local
         % variable; these get added to the map so that proper constraints
         % can be setup later
         line = size (new_edges, 1);
         map = [map; new_point, line];
      end;

      % create a line from the last (possible interior point) to the end of
      % the line, completing the number of segments. if there were no
      % interior points, then this will represent the line from the old set
      new_edges = [new_edges; last_point, old_edges(i, 2), c];
   end;