function G = computeVEM2DGeometry(G)
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

% x = G.nodes.coords(G.cells.nodes,1);
% y = G.nodes.coords(G.cells.nodes,2);
% ncn = diff(G.cells.nodePos);
% [ii,jj] = blockDiagIndex(ncn);
% xx = sparse(ii,jj,rldecode(x,rldecode(ncn,ncn,1),1));
% yy = sparse(ii,jj,rldecode(y,rldecode(ncn,ncn,1),1));
% bsxfun(@minus,xx,x)
% bsxfun(@minus,yy,y)

for i = 1:G.cells.num
    nodeNum = G.cells.nodePos(i):G.cells.nodePos(i+1)-1;
    nodes = G.cells.nodes(nodeNum);
    X = G.nodes.coords(nodes,:);
    cellDiameters(i) = cellDiameter(X);
end

G.cells.('diameters') = cellDiameters;

end