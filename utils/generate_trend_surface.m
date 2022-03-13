function [x, y, zgrid, surf] = generate_trend_surface(points, num_coefs, ...
                                                     samples, stiffness)
%
% Generate a spline approximating surface ("trend surface") from a set of
% scattered points, and sample it into a regular grid.
%
% SYNOPSIS:
%   function [x, y, zgrid, surf] = generate_trend_surface(points, num_coefs, samples, stiffness)
%
% PARAMETERS:
%   points    - Nx3 array with coordinates of N scattered 3D-points (must be
%               projectable to a surface in the (x, y) plane).
%   num_coefs - Number of knot intervals for the x-cooridinate (u-parameter)
%   samples   - resolution (x and y direction) of sampled regular grid
%   stiffness - stiffness parameter, defining the tradeoff between surface
%               regularity and approximation quality
% RETURNS:
%   x     - x-coordinates for zgrid
%   y     - y-coordinates for zgrid
%   zgrid - sampled grid from approximating surface
%   surf  - spline surface approximating the scattered points

%{
Copyright 2009-2022 SINTEF Digital, Mathematics & Cybernetics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}

   % fitting trend surface
   surf = BivariateSplineFunction([4, 4], num_coefs, num_coefs);
   surf.approximate(points(:,1:2), points(:,3), 'data', stiffness);

   % sampling trendsurface
   box = [min(points(:,1)), min(points(:,2)), max(points(:,1)), max(points(:,2))];
   
   x = linspace(box(1), box(3), samples);
   y = linspace(box(2), box(4), samples);
   
   [U, V] = ndgrid(x, y);

   zgrid = reshape(surf.evaluate([U(:), V(:)], 0, 0), samples, samples);
end
