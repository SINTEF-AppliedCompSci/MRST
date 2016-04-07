G = cartGrid([4,4]);

G = computeGeometry(G);
G = mrstGridWithFullMappings(G);

b = 1:G.faces.num;
b = b(G.faces.centroids(:,1) == 2)

VEM2D_makeInternalBoundary(G,b);