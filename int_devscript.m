mrstModule add vem vemmech

%%

n = 2;
G = computeGeometry(cartGrid([n,n], [2,2]));
G.nodes.coords = G.nodes.coords;
G = computeVEMGeometry(G);

p = Polynomial([2,0; 1,1], [0.4; 1]);
val = integratePolygon(p, G, (1:G.cells.num)', 6);

%%

degree = 2;
C = computeCellIntegrals(G, degree)