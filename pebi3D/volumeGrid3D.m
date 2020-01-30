function [grid_3]  = volumeGrid3D(pdim, faults, grids_2, ds, gamma)
% Find the sites and 3D grid that conform to surfaces, and the intersection
% of surfaces.
%
% SYNOPSIS:
%   grids_3 = volumeGrid3D(pdim, fautls, grids_2, ds, gamma)
%
% PARAMETERS
%   pdim            - Vector, length numel(celldim), of physical size in
%                     units of meters of the computational domain.
%   fautls          - cell array of faults. Each element is a vector of
%                     the vertices of the fault
%   grids_2         - cell array of 2D grids as returned from
%                     surfaceGrid3D
%   ds              - Cell size of the 2D cells
%   gamma           - Distance surface sites are placed from the 1D grids
%
% RETURNS:
%   grids_3         - The 3D conforming MRST-grid. In
%                     addition to a normal MRST-grid the field
%                     G.cells.sites is added, which is an array of the
%                     sites used to generate the grid
%
% EXAMPLE:
%   f1 = [1,3,2; 4,3,2; 4,3,4; 1,3, 4];
%   f2 = [2,2,3.3; 5,2,3.3; 5,4,3.3; 2,4, 3.3];
%   fracs = {f1, f2};
%   intersections = surfaceIntersections3D(fracs);
%   grids_1 = lineGrid3D(intersections, 0.2);
%   grids_2 = fautlSites(fracs, grids_1, intersections, 0.2, 0.2/6)
%   grids_3 = volumeGrid3D([6,6,6], fracs, grids_2, 0.4, 0.1)
%   plotGrid(grids_3)
%
% SEE ALSO
%   lineGrid3D, surfaceIntersections3D, surfaceGrid3D, compositePebiGrid2D, pebi, surfaceSites2D, lineSites2D.

%{
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2016 Runar Lie Berge. See COPYRIGHT.TXT for details.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%}  
internal_pts = zeros(0, 3);

for f = 1:numel(grids_2)
    f_pts = faults{f};
    normal = normalFromPoints(f_pts);
    pts_2 = grids_2{f}.cells.sites; 

    pts_3 = bsxfun(@plus, pts_2, normal * gamma);
    pts_3_neg = bsxfun(@minus, pts_2, normal * gamma);
    internal_pts = [internal_pts; pts_3; pts_3_neg];
end


% Create reservoir sites
x = pdim(1); y = pdim(2); z = pdim(3);

ds = ds * 1.2;
nx = x / ds;
ny = y / ds;
nz = z / ds;

xa = linspace(ds/2, x-ds/2, nx + 1);
ya = linspace(ds/2, y-ds/2, ny + 1);
za = linspace(ds/2, z-ds/2, nz + 1);

bdr   = [ 0, 0, 0;  ...
          x, 0, 0;  ...
          x, y, 0;  ...
          0, y, 0;  ...
          0, 0, z;  ...
          x, 0, z;  ...
          x, y, z;  ...
          0, y, z];

[X,Y,Z] = ndgrid(xa,ya,za);
rSites = [X(:), Y(:), Z(:)];

sites = cellfun(@(c) c.cells.sites, grids_2, 'un', false);
CC = vertcat(sites{:});

rSites = surfaceSufCondFromGrid3D(rSites, grids_2, gamma);
sites_3 = [internal_pts; rSites];

grid_3 = mirroredPebi3D(sites_3, bdr);
grid_3.cells.sites = sites_3;


end
