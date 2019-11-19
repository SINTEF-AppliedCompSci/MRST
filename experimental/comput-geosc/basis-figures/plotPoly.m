mrstModule add vem vemmech upr

%%

pth = fullfile(mrstPath('dg'), 'examples', 'comput-geosc', 'basis-figures', 'fig');
saveeps = @(name) print(fullfile(pth, name), '-depsc');

%%

close all

rng(1);
G = computeVEMGeometry(pebiGrid(1/5, [1,1]));
G.cells.order = (1:G.cells.num)';
plotToolbar(G, G);

%%

R = @(t)[cos(t), sin(t); -sin(t), cos(t)];

figure
cellNo = 22;
n = G.cells.nodes(G.cells.nodePos(cellNo):G.cells.nodePos(cellNo+1)-1);
x = G.nodes.coords(n,:);

x = x*R(pi/8);

gr = [1,1,1]*0.8;
patch(x(:,1), x(:,2), gr);

axis equal tight
box on
ax = gca;
[ax.XTick, ax.YTick] = deal([]);

saveeps('cell-2d')

%%

close all


xmax = [1,1,1];
%     xmax = [500, 1000, 500]/2*meter;
%     xmax = [100, 300, 100]*meter;
%     n = [5, 10, 5];
n = [5,5,5];

xx = cell(3,1);
d = 5*meter/1000;
for dNo = 1:3
    xx{dNo} = linspace(d,xmax(dNo)-d,n(dNo));
end

[x,y,z] = ndgrid(xx{:});
x = [x(:), y(:), z(:)];
npts = size(x,1);

rng(1);
d = xmax./n*0.25;
x = x + randn(npts,3).*d;

bnd = [0      , 0      , 0      ;
       xmax(1), 0      , 0      ;
       xmax(1), xmax(2), 0      ;
       0      , xmax(2), 0      ;
       0      , 0      , xmax(3);
       xmax(1), 0      , xmax(3);
       xmax(1), xmax(2), xmax(3);
       0      , xmax(2), xmax(3)];
bnd = delaunayTriangulation(bnd);
G = clippedPebi3D(x, bnd);

G = computeVEMGeometry(G);
% G = computeCellDimensions(G);

plotGrid(G)

%%

R = @(t)[cos(t), sin(t); -sin(t), cos(t)];

figure
[f, c] = boundaryFaces(G);
isInt = true(G.cells.num,1);
isInt(c) = false;
intCells = find(isInt);
cellNo = intCells(3);

% x = x*R(pi/8);

gr = [1,1,1]*0.8;
plotGrid(G, cellNo, 'facec', gr);

axis equal tight
box on
ax = gca;
[ax.XTick, ax.YTick, ax.ZTick] = deal([]);
view([127,21])

saveeps('cell-3d')
