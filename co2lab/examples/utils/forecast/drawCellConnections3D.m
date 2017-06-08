function drawCellConnections3D(Gt, cells, varargin)

    x = Gt.cells.centroids(cells, 1);
    y = Gt.cells.centroids(cells, 2);
    z = Gt.cells.z(cells);
    smooth=@(x) x;
    x = smooth(x);
    y = smooth(y);
    z = smooth(z);
    plot3(x, y, z, varargin{:});
    
end
