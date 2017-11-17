function [grid_3]  = reservoir_sites(pdim, faults, grids_2, ds, gamma)

internal_pts = zeros(0, 3);

for f = 1:numel(grids_2)
    f_pts = faults{f};
    normal = normal_from_points(f_pts);
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
R = ds * ones(size(CC, 1), 1);

%rSites = faultSufCond3D(rSites, CC, R);
rSites = faultSufCondFromGrid3D(rSites, grids_2, gamma);
sites_3 = [internal_pts; rSites];
%grid_3 = CVD3D(rSites, bdr, 'fixedPts', internal_pts, 'maxIt', 30);

grid_3 = mirroredPebi(sites_3, bdr);
grid_3.cells.sites = sites_3;
%grid_3 = internal_pts;

end