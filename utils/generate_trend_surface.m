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
%
% EXAMPLE:
%
% SEE ALSO:
%

   

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