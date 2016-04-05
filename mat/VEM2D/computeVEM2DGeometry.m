function G = computeVEM2DGeometry(G,f,k, alpha)

G = computeGeometry(G);
G = mrstGridWithFullMappings(G);
G = sortEdges(G);

cellDiameters = zeros(G.cells.num,1);
for i = 1:G.cells.num
    nodeNum = G.cells.nodePos(i):G.cells.nodePos(i+1)-1;
    nodes = G.cells.nodes(nodeNum);
    X = G.nodes.coords(nodes,:);
    cellDiameters(i) = cellDiameter(X);
end

G.cells.('diameters') = cellDiameters;

% [AK, bK, SK, PN] = VEM2D_projectors(G, f, k, alpha);
% 
% G.cells.('AK') = AK;
% G.cells.('bK') = bK;
% G.cells.('SK') = SK;
% G.cells.('PN') = PN;


end