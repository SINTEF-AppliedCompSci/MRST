function [monomialNodeVals, monomialEdgeVals] = monomialValues(G)

m3D =      @(X) [ones(size(X,1),1) , ...
                X(:,1)              , ...   %   (1,0,0)
               X(:,2)               , ...   %   (0,1,0)
               X(:,3)               , ...   %   (0,0,1)
               X(:,1).^2          , ...   %   (2,0,0)
               X(:,1).*X(:,2) , ...   %   (1,1,0)
               X(:,1).*X(:,3), ...   %   (1,0,1)
               X(:,2).^2, ...   %   (0,2,0) 
               X(:,2).*X(:,3), ...   %   (0,1,1)
               X(:,3).^2];      %   (0,0,2)icenter of K.

nodeNum = mcolon(G.cells.nodePos(1:end-1),G.cells.nodePos(2:end)-1);
nodes = G.cells.nodes(nodeNum);
X = G.nodes.coords(nodes,:);

edgeNum = mcolon(G.cells.edgePos(1:end-1),G.cells.edgePos(2:end)-1);
edges = G.cells.edges(edgeNum);
Ec = G.edges.centroids(edges,:);

X = bsxfun(@rdivide, X - rldecode(G.cells.centroids, diff(G.cells.nodePos), 1), ...
                         rldecode(G.cells.diameters, diff(G.cells.nodePos), 1));
Ec = bsxfun(@rdivide, Ec- rldecode(G.cells.centroids, diff(G.cells.edgePos), 1), ...
                         rldecode(G.cells.diameters, diff(G.cells.edgePos), 1));

monomialNodeVals = m3D(X);
monomialEdgeVals = m3D(Ec);

end