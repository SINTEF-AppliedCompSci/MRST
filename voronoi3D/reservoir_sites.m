function grid_3  = reservoir_sites(pdim, faults, grids_2, ds)

internal_pts = zeros(0, 3);
gamma = 1/8;
for f = 1:size(grids_2, 1)
    f_pts = faults{grids_2{f, 2}};
    normal = normal_from_points(f_pts);
    pts_2 = grids_2{f, 1};

    pts_3 = bsxfun(@plus, pts_2, normal * ds * gamma);
    pts_3_neg = bsxfun(@minus, pts_2, normal * ds * gamma);
    internal_pts = [internal_pts; pts_3; pts_3_neg];
end


% Create reservoir sites
x = pdim(1); y = pdim(2); z = pdim(3);

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

CC = vertcat(grids_2{:,1});
R = ds * ones(size(CC, 1), 1);

rSites = faultSufCond3D(rSites, CC, R);

pts = [internal_pts; rSites];

grid_3 = mirroredPebi(pts, bdr);
%grid_3 = internal_pts;

end