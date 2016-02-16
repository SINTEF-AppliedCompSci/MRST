function G = tetrahedronCube(dims, gridLim, ran)

nPx = dims(1)+1;
nPy = dims(2)+1;
nPz = dims(3)+1;

dx = gridLim(1)/(nPx-1);
xx = 0:dx:gridLim(1);
dy = gridLim(2)/(nPy-1);
yy = 0:dy:gridLim(2);
dz = gridLim(3)/(nPz-1);
zz = 0:dz:gridLim(3);

[x, y, z] = meshgrid(xx, yy, zz);   

if ran == 1
    assert(min(dims) > 1, ...
           'Random grid requires 2 or more cells in each coordinate direction.')
    x(:,2:nPx-1,:) = x(:,2:nPx-1,:) + random('Normal', 0, dx/4, nPy, nPx-2, nPz);
    y(2:nPy-1,:,:) = y(2:nPy-1,:,:) + random('Normal', 0, dy/4, nPy-2, nPx, nPz);
    z(:,:,2:nPz-1) = z(:,:,2:nPz-1) + random('Normal', 0, dz/4, nPy, nPx, nPz-2);
end

P = [x(:), y(:), z(:)];
T = DelaunayTri(P);

G = tetrahedralGrid(P, T.Triangulation);

end