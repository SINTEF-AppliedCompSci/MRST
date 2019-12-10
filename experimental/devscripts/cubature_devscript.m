n = 10;
G = computeGeometry(cartGrid([n,n], [10,10]));

cubTri = TriangleCubature(G, 2);
cubLin = LineCubature(G, 2);

n = 3;
G = computeVEMGeometry(cartGrid([n,n,n], [10,10,10]));
cubTet = TetrahedronCubature(G, 2);

close all
figure
hold on
plotGrid(G, 'facec', 'none')
x = cubTet.points;
plot3(x(:,1), x(:,2), x(:,3), '.');

