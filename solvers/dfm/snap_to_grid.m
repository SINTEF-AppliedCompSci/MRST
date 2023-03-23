function [vertices] = snap_to_grid(vertices, box, opts)
%
% Move vertices to the closest (structured) grid point.
%
% Use this function to avoid having more than one point within each (fine)
% grid cell. If this is not done, a low quality grid may result
%
% Copyright (C) 2006-2007 Uni Research AS
%
% This file is licensed under the GNU General Public License v3.0.
%

  x = 1; y = 2; minimum = 1; maximum = 2;

  % Precision is the grid size of the underlying structured grid
  precision = opts.precision;

  if precision > eps && ~isempty (vertices),
     horz_prec = (box(maximum, x) - box(minimum, x)) .* precision;
     vert_prec = (box(maximum, y) - box(minimum, y)) .* precision;
     vertices = [round(vertices (:, x) ./ horz_prec) .* horz_prec, ...
         round(vertices (:, y) ./ vert_prec) .* vert_prec];
  end;