function [pt] = lines_intersect(x_a_1,y_a_1,x_b_1,y_b_1,x_a_2,y_a_2,x_b_2,y_b_2)
% [pt] = lines_intersect (vertices, edge1, edge2)
%
% Determine the intersection between two lines, if any. t will either be
% 0 or 1 depending on whether there are any intersections, in that case x
% and y will be the coordinates of the point that intersect. The input
% arguments is the global list of points, and two tuples that contains
% the indices into this list.
%
% Portions Copyright (C) 2006-2007 Uni Research AS
% This file is licensed under the GNU General Public License v3.0.

  % dissect the input into individual coordinates. index into edgeX array
  % is the point (a or b) on line X. second index into vertices is the
  % coordinate (x or y) for that point. read the names of these variables
  % as "x for the point a on line 1". (The use of the return variables
  % not intended; they are used as named constants here).

  % these should have some values, although they are unspecified if no
  % match is found in the later tests
  pt = [];

  % bounding boxen for each of the two segments
  x_min_1 = min (x_a_1, x_b_1);
  x_max_1 = max (x_a_1, x_b_1);
  y_min_1 = min (y_a_1, y_b_1);
  y_max_1 = max (y_a_1, y_b_1);
  x_min_2 = min (x_a_2, x_b_2);
  x_max_2 = max (x_a_2, x_b_2);
  y_min_2 = min (y_a_2, y_b_2);
  y_max_2 = max (y_a_2, y_b_2);

  % the two lines may only intersect if their bounding boxen
  % intersect. we don't mind if their extensions intersect. test each
  % dimension the same way (one line per dimension) and then repeat the
  % test for the case where the segments have opposite roles.
  proximity = ( ( x_min_2 <= x_max_1 && x_max_2 >= x_min_1 ) && ...
                ( y_min_2 <= y_max_1 && y_max_2 >= y_min_1 ) ) || ...
              ( ( x_min_1 <= x_max_2 && x_max_1 >= x_min_2 ) && ...
                ( y_min_1 <= y_max_2 && y_max_1 >= y_min_2 ) );
  if proximity
    % now we know that they are in the proximity of eachother; check if
    % they actually cross.

    % check on which side (-1, 0 or 1) that the first point of the second
    % line is in relation to the entire first line, by calculating the
    % cross-product of the vector from the free point to the line and the
    % vector of the line itself.
    det_1_a = (x_b_1 - x_a_1) * (y_a_2 - y_a_1) - ...
              (x_a_2 - x_a_1) * (y_b_1 - y_a_1);

    % same thing, but with the second point of the second line
    det_1_b = (x_b_1 - x_a_1) * (y_b_2 - y_a_1) - ...
              (x_b_2 - x_a_1) * (y_b_1 - y_a_1);

    % symmetric, but with points from line 1 measured against line 2
    det_2_a = (x_b_2 - x_a_2) * (y_a_1 - y_a_2) - ...
              (x_a_1 - x_a_2) * (y_b_2 - y_a_2);

    det_2_b = (x_b_2 - x_a_2) * (y_b_1 - y_a_2) - ...
              (x_b_1 - x_a_2) * (y_b_2 - y_a_2);

    % if these two have different signs, then the two lines cross each
    % other (if one line crosses the other, then the other must also cross
    % that one!), or if one of the points actually is on the other line
    straddle = ( sign (det_1_a * det_1_b) <= 0 && ...
                 sign (det_2_a * det_2_b) <= 0 );

    if straddle
      % represent each line on the form p*x + q*y = r
      p_1 = y_a_1 - y_b_1;
      q_1 = x_b_1 - x_a_1;
      r_1 = p_1 * x_a_1 + q_1 * y_a_1;

      % same as above, but for the second line
      p_2 = y_a_2 - y_b_2;
      q_2 = x_b_2 - x_a_2;
      r_2 = p_2 * x_a_2 + q_2 * y_a_2;

      % multiply the equation for the first line by the second
      % coefficient (q) of the second line and vice versa (second
      % equation is also negated to get the same denominator). subtract the
      % equations and solve for each variable.
      det = p_1 * q_2 - p_2 * q_1;
      if abs (det) > 1e-6

        x = (r_1 * q_2 - r_2 * q_1) ./ det;
        y = (p_1 * r_2 - p_2 * r_1) ./ det;

        % here we could have decided that we didn't want to return
        % endpoints as intersections, but there are too many special cases
        % to think about, so it is more efficient to simple leave those
        % problems to the split procedure (it gets to decide).

        % return a full tuple that tells us that they cross
        pt = [x y];
      end;
    end;
  end;
