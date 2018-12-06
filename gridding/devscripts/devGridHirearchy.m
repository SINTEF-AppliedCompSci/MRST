mrstModule add vemmech coarsegrid msrsb

%%

n = 20;
G = pebiGrid(1/n, [1,1]);

plotGrid(G)
axis equal tight

%%

G = computeGeometry(G);
G = computeCellDimensions2(G);

%%

% Cartesian coarse grid
G_cart = cartGrid([100, 100]);
p_cart = partitionUI(G_cart, [10, 10]);
p_cart = sampleFromBox(G, reshape(p_cart, G_cart.cartDims));
GC = generateCoarseGrid(G, p_cart);
GC = coarsenGeometry(GC);
GC = addCoarseCenterPoints(GC);

%%

GC = coarsenCellDimensions(GC);

%%

close all

cNo = 45;
plotGrid(GC)
plotGrid(GC, cNo, 'facec', 'r');
hold on
x = [GC.cells.xMin(cNo,:);GC.cells.xMax(cNo,:)];
plot([x(1,1), x(2,1), x(2,1), x(1,1), x(1,1)], ...
     [x(1,2), x(1,2), x(2,2), x(2,2), x(1,2)], 'b', 'linew', 2);
 
f = GC.cells.faces(GC.cells.facePos(cNo):GC.cells.facePos(cNo+1)-1);
for fNo = 1:numel(f)
    x = GC.nodes.coords(GC.faces.nodePos(f(fNo)):GC.faces.nodePos(f(fNo)+1)-1,:);
    plot(x(:,1), x(:,2), '--g', 'linew', 2);
end
axis equal tight

%%

n = 5;
G = computeGeometry(cartGrid([n,n,n], [1,1,1]));
t = pi/8;
R = [cos(t), -sin(t), 0; sin(t), cos(t), 0; 0, 0, 1];
G.nodes.coords = G.nodes.coords*R';

close all
plotGrid(G)
axis equal tight

%%

G = computeGeometry(G);
G = computeCellDimensions2(G);

%%

close all

cNo = 11;
plotGrid(G, 'facec', 'none')
plotGrid(G, cNo, 'facec', 'r');
hold on
x = [G.cells.xMin(cNo,:);G.cells.xMax(cNo,:)];

G.cells.dx(cNo,:)

faces = G.cells.faces(G.cells.facePos(cNo):G.cells.facePos(cNo+1)-1);
G.faces.dx

view(3)
axis equal tight