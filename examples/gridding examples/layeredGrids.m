%% Example of creating layered grids
%
% This example demonstrates the use of 'makeLayeredGridNWM' to extrude 2D 
% grid to the layered 3D grid according the topology of 2D grid and 
% provided layer point set. The layer points are given on all surfaces, and
% topologically aligned in layered direction.

mrstModule add nwm

%% Make a layered Cartesian grid
% Build the 2D grid
g = cartGrid([20, 20], [200, 200]);
% The layer point set
xy = g.nodes.coords;
n = size(xy, 1);
Z = linspace(0, 10, 10);
players = arrayfun(@(z)[xy, z*ones(n,1)], Z, 'UniformOutput', false);
% Make the layered grid
G = makeLayeredGridNWM(g, players);

% Plot the grid. We can plot the specified layers / surfaces by calling 
% 'G.cells.layers==layer'/ 'G.faces.surfaces==surface'
figure, hold on; axis off
plotGrid(G, 'facecolor', 'none')
plotGrid(G, G.cells.layers == 3, 'facecolor', [.0, .44, .74])
plotFaces(G, G.faces.surfaces == 8, 'facecolor',  [.85, .32, .09])
legend('Layered Cartesian grid', 'Cells on layer 3', 'Faces on surface 8')
view(3)

%% Make a layered radial-Cartesian hybrid grid with inclination
% Build the Cartesian grid
GC = cartGrid([15, 15], [300, 300]);
GC = computeGeometry(GC);
% Define the well region by logical indices
ij = gridLogicalIndices(GC);
idxI = ij{1} >= 6 & ij{1} <= 10 & ij{2} >= 6 & ij{2} <= 10;
cI = find(idxI); 
% Place the well at the region center
pW  = [150, 150];
% Define radial parameters
[nR, rW, rM] = deal(10, 0.2, 35);
% Get the hybrid grid
g = radCartHybridGrid(GC, cI, rW, rM, nR, pW);

% The layer point set, let the surfaces inclines in X direction
x  = g.nodes.coords(:,1);
dx = linspace(0, 100, 12);
y = g.nodes.coords(:,2);
n = size(x, 1);
Z = linspace(10, 100, 12);
players = arrayfun(@(dx, z)[x-dx, y, z*ones(n,1)], dx, Z, 'UniformOutput', false);
% Make the layered grid
G = makeLayeredGridNWM(g, players);

% Plot the inclined grid
figure, hold on; axis off
plotGrid(G, 'facecolor', [.0, .44, .74])
view(3)

% Plot the radial subgrid
cells = reshape((1:G.cells.num)', g.cells.num, G.layers.num);
radDims = g.subGrids{1}.radDims;
radIdx = ( 1: (radDims(1)*(radDims(2)-1)) )';
cellsRad = cells(radIdx, :);
figure, hold on; axis off
plotGrid(G, 'facecolor', 'none')
plotGrid(G, cellsRad, 'facecolor', [.85, .32, .09])
view(3)

%% Make a layered horizontal radial grid
[nA, nR, rW, rM] = deal(40, 10, 1, 10);
th = linspace(0, 2*pi, nA+1); th = th(1:end-1);
r = logspace(log10(rW), log10(rM), nR+1);
[R, TH] = meshgrid(r, th);
[px, py] = pol2cart(TH(:), R(:));
p = [px(:), py(:)];
[g, t] = buildRadialGrid(p, nA, nR);

% The layer point set, convert the grid from XY plane to YZ plane
yz = g.nodes.coords;
n = size(yz, 1);
X = linspace(0, 200, 12);
players = arrayfun(@(x)[x*ones(n,1), yz], X, 'UniformOutput', false);
% Make the layered grid. Note if the 2D grid are not on the XY plane, the
% connectivity list of the 2D grid should be provided.
G = makeLayeredGridNWM(g, players, 'connectivity', t);

% Plot the grid. 
figure, hold on; axis equal tight off
plotGrid(G, 'facecolor', 'none')
plotGrid(G, G.cells.layers == 3, 'facecolor', [.0, .44, .74])
plotFaces(G, G.faces.surfaces == 8, 'facecolor',  [.85, .32, .09])
legend('Layered Cartesian grid', 'Cells on layer 3', 'Faces on surface 8')
view(3)

