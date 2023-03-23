function grids_1 = lineGrid3D(intersections, ds)
% Find the sites of the grids defined by intersection of surfaces, and generate
% the corresponding 1D grids
%
% SYNOPSIS:
%   grids_1 = lineGrid3D(intersections, ds)
%
% PARAMETERS
%   intersecitons   - The intersection of surfaces as returned from
%                     surfaceIntersections3D
%   ds              - Cell size of the 1D cells
%
% RETURNS:
%   grids_1         - cell-array where each element is a 1D grid. The grids
%                     are defined by a set of vertices and the
%                     corresponding PEBI-site
%
% EXAMPLE:
%   f1 = [1,3,2; 4,3,2; 4,3,4; 1,3, 4];
%   f2 = [2,2,3.3; 5,2,3.3; 5,4,3.3; 2,4, 3.3];
%   fracs = {f1, f2};
%   intersections = surfaceIntersections3D(fracs);
%   grids_1 = lineGrid3D(intersections, 0.2);
%   % the PEBI-sites are:
%   grids_1{1}{1} 
%   % the vertices of the grid are:
%   grids_1{1}{2}
% SEE ALSO
%   lineGrid3D, surfaceGrid3D, volumeGrid3D, compositePebiGrid2D, pebi, surfaceSites2D, lineSites2D.

%{
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2016 Runar Lie Berge. See COPYRIGHT.TXT for details.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%}  
grids_1 = {};
for i = 1:numel(intersections)
   int = intersections{i} ;
   p1 = int{1}(1,:);
   p2 = int{1}(2,:);
   v = p2 - p1;
   L = sqrt(sum(v.^2));
   v = v / L;
   n = round(L / ds);
   if L < 1.5*ds
       pts = zeros(0, 3);
       vert = zeros(0, 3);
   else
       t = linspace(ds/2, L - ds/2, n);
       t_vert = [0, (t(1:end-1) + t(2:end))/2, L];
       steps = bsxfun(@times, v, t');
       steps_vert = bsxfun(@times, v, t_vert');
       pts = bsxfun(@plus, p1, steps);
       vert = bsxfun(@plus, p1, steps_vert);
   end
   grids_1 = {grids_1{:}, {pts, vert}};
end


end
