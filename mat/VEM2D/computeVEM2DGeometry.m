function G = computeVEM2DGeometry(G)
%--------------------------------------------------------------------------
%   Computes VEM geometry of MRST grid G.
%
%   SYNOPSIS:
%       G = computeVEM2DGeometry(G)
%
%   DESCRIPTION:
%       Computes geometry using MRST functions G = computeGeometry(G) and
%       G = mrstGridWithFullMappings(G), and computes cell diameters.
%       See MRTS functions computeGeomerty and mrstGridWithFullMappings for
%       details and copyright info.
%
%   REQUIRED PARAMETERS:
%       G   - 2D MRST grid.
%
%   RETURNS:
%       G   - Grid with computed VEM geometry.  
%-----------------------------------------------------------------ØSK-2016-

assert(G.griddim == 2, 'Function only supprots 2D grids');

G = computeGeometry(G);
G = mrstGridWithFullMappings(G);

cellDiameters = zeros(G.cells.num,1);
for i = 1:G.cells.num
    nodeNum = G.cells.nodePos(i):G.cells.nodePos(i+1)-1;
    nodes = G.cells.nodes(nodeNum);
    X = G.nodes.coords(nodes,:);
    cellDiameters(i) = cellDiameter(X);
end

G.cells.('diameters') = cellDiameters;

end