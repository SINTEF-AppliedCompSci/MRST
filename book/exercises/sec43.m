%% Extend triangleGrid to triangulated surfaces
load trimesh3d;
G = myTriangleGrid([x y z], tri); clf, plotGrid(G); view(3);

G = computeGeometry(G);
hold on,
c = G.cells.centroids;
plot3(c(:,1),c(:,2),c(:,3),'.r');
f = G.faces.centroids;
plot3(c(:,1),c(:,2),c(:,3),'.b');
hold off
