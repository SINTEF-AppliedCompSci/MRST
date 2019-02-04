%% Exercise 3.4.1
% Revisit Exercise 3.1.1, but use lower grid resolution 
dx = .08;
x = 1-0.5*cos((-1:dx:1)*pi);
x = -1 - 1.5*dx+dx*cumsum(x);
y = 0:0.1:.5;
z = 0:0.08:1;
G = computeGeometry(tensorGrid(x, y, z));
cyl = @(x) x(:,1).^2 + x(:,3).^2;
G = removeCells(G,cyl(bsxfun(@minus, G.cells.centroids,[0 0 0]))<0.16);
G = removeCells(G,cyl(bsxfun(@minus, G.cells.centroids,[0 0 1]))<0.16);
clf, subplot(1,2,1)
plotCellData(G,G.cells.volumes,'EdgeAlpha',.2); view(-20,25); axis equal tight
cx = caxis;

% Perturb all points except for those on the exerior surface and also
% have a normal vector with no component in the y-direction
ofaces = boundaryFaces(G);
zeroy  = abs(G.faces.normals(ofaces,2))>eps;
pflag  = true(G.nodes.num,1);
npert  = mcolon(G.faces.nodePos(ofaces(zeroy)), G.faces.nodePos(ofaces(zeroy)+1)-1);
npert  = unique(G.faces.nodes(npert));
pflag(npert) = false;
pts = G.nodes.coords;
pts = pts + bsxfun(@times, randn(G.nodes.num,3)*.01, double(pflag));
G.nodes.coords = pts;
G = computeGeometry(G);
subplot(1,2,2)
plotCellData(G,G.cells.volumes,'EdgeAlpha',.2); view(-20,25); axis equal tight
caxis(cx);

%% Exercise 3.4.2
% Simple check to see that the extended triangleGrid function from Exercise
% 3.2.2 gives reasonable cell and face centroids. First in 2D
load trimesh2d;
G = computeGeometry(myTriangleGrid([xfe 0*xfe yfe], trife));
g = computeGeometry(triangleGrid([xfe yfe],trife));

clf
subplot(2,2,1), 
plot(G.cells.volumes, g.cells.volumes,'o');
v = [G.cells.volumes; g.cells.volumes];
hold on
plot([min(v) max(v)], [min(v) max(v)],'r-');
hold off
axis tight

subplot(2,2,3), plot(G.faces.areas, g.faces.areas, 'o');
a = [G.faces.areas; g.faces.areas];
hold on
plot([min(a) max(a)], [min(a) max(a)],'r-');
hold off, axis tight

subplot(2,2,[2 4]),
plotGrid(G); view(3);
hold on,
c = G.cells.centroids;
plot3(c(:,1),c(:,2),c(:,3),'.r');
f = G.faces.centroids;
plot3(f(:,1),f(:,2),f(:,3),'.b');
hold off
view(0,0); axis([0 100 0 1 0 100]);

%%
% Then we do a similar test in 3D
load trimesh3d;
G = myTriangleGrid([x y z], tri);
G = computeGeometry(G);

clf, plotGrid(G,'Facealpha',1); view(3); axis tight
hold on,
c = G.cells.centroids;
plot3(c(:,1),c(:,2),c(:,3),'.r');
f = G.faces.centroids;
plot3(f(:,1),f(:,2),f(:,3),'.b');
hold off

%% Exercise 3.4.3
% Create grid with invalid nodes
G = cartGrid([10 10]);
G.nodes.coords([10 15 20],1)=NaN;
G.nodes.coords(30+[10 15 20],2)=NaN;
clf, subplot(1,2,1)
plotGrid(G);
title(['# cells: ' num2str(G.cells.num)]);

% Purge invalid grid cells
cn = cellNodes(G);
nflag = ~any(isnan(G.nodes.coords),2);
cflag = accumarray(cn(:,1), nflag(cn(:,3))) < accumarray(cn(:,1),1);
g = extractSubgrid(G,~cflag);
% Alternative: g = removeCells(G,cflag);
subplot(1,2,2)
plotGrid(g);
title(['# cells: ' num2str(g.cells.num)]);
