function [G, cornercelltbl] = gridForPaperConvTest(Nx, gridType)
% We split the square (cube in 3d) in 4 (8 in 3d) parts. Indices of the cells of
% the (upper) north-west are returned.
    
% Planned Grid type:
% 1: Cartesian
% 2: Triangles by alternating bisection of triangles
% 3: Equilateral triangles
% 4: Triangles by uniform bisection (grid greated by Dolfin)
    
    Nd = numel(Nx);
    
    switch gridType
      case 1
        G = computeGeometry(cartGrid(Nx, ones(1, Nd)));
        c = G.cells.centroids;
        cornercelltbl.cells = find(all(c > 0.5, 2));
        cornercelltbl = IndexTable(cornercelltbl);
      otherwise
        error('gridType not recognized');
    end
end
