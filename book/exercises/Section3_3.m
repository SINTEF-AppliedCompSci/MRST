%% Exercise 3.3.1
% 
dx = .04;
x = 1-0.5*cos((-1:dx:1)*pi);
x = -1 - 1.5*dx+dx*cumsum(x);
y = 0:0.05:.5;
z = 0:0.04:1;
G = computeGeometry(tensorGrid(x, y, z));
cyl = @(x) x(:,1).^2 + x(:,3).^2;
G = removeCells(G,cyl(bsxfun(@minus, G.cells.centroids,[0 0 0]))<0.16);
G = removeCells(G,cyl(bsxfun(@minus, G.cells.centroids,[0 0 1]))<0.16);
clf, plotCellData(G,G.cells.volumes,'EdgeAlpha',.2); view(-20,25); axis equal tight


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
