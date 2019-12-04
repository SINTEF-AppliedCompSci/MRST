%% Demonstration of building an individual radial grid
%
% The radial grid offers superior performance in describing the radial flow
% in the well vicinity. The mrst grid factory has provided the grid
% constructor 'pebi' which can build the radial grid from a set of radially
% distributed points (see 'showPEBI' in 'modules\book\examples\grids'). 
% However, this method has several restrictions. In this example, we show
% the restrictions of the pebi-style radial grid, and then demonstrate how 
% to build an radial grid using the constructor 'tessellationGrid'. The 
% generating points are the grid nodes, and the indices of grid objects, 
% i.e. cells, faces, and nodes, are numbered logically.

clear
mrstModule add nwm

%% Generate radial point set
% nA - Number of cells in angular direction
% nR - Number of cells in radial direction
% rW - The minimum radius
% rM - The maximum radius

[nA, nR, rW, rM] = deal(40, 10, 1, 10);
th = linspace(0, 2*pi, nA+1); th = th(1:end-1);
r = logspace(log10(rW), log10(rM), nR+1);
[R, TH] = meshgrid(r, th);
[px, py] = pol2cart(TH(:), R(:));
% The points obey the logical numbering (angularly cycle fastest, then 
% radially)
p = [px(:), py(:)];

%% Build the radial grid by 'pebi'
% Note 'pebi' requires the continuous region, i.e. no 'hole' inside the
% region. So, the central points should be added into the generation point
% set.
pP = [[0, 0]; p];
H = pebi(triangleGrid(pP));
H = computeGeometry(H);

% Remove the central cell
d = sqrt( sum(H.cells.centroids.^2, 2) );
cC = find(d == min(d));
H = removeCells(H, cC);

% 1. The generating points will not be the grid nodes, instead they are
% located inside the grid cells (see the definition of the Voronoi
% diagram), which goes against the situations where we want to locate the 
% specific points, e.g. gluing the grids, modeling the wellbore explicitly.
figure, hold on; axis equal off
plotGrid(H)
demoPlotLine(p, 'ks', 'r', 2)
title('Pebi-style radial grid and generating points')

% 2. The pebi style radial grid is still an unstructured grid. Some cells 
% have more than 4 faces. For another, conflict points (point too close to 
% each other) always appear in such grids, resulting in very small faces. 
% Some postprocessing, e.g. 'removeShortEdges', is required.

% Display the face numbers of cells
cfN = diff(H.cells.facePos);
figure, hold on; axis equal off
plotCellData(H, cfN), colorbar
title('Face numbers of cells in Pebi-style radial grid')

% Display the number of very small faces
fprintf(' There are altogether %d faces with areas smaller than 1e-10\n', ...
    nnz( H.faces.areas < 1e-10 ))

% 3. The indices of faces are not numbered logically. We cannot find the
% radial faces and angular faces directly. This would be unfavorable when
% computing the radial transmissibility. 
c = find(cfN==4);
c = c([5:10:end]);
figure, hold on; axis equal off
plotGrid(H, 'facecolor', 'none')
cols = {'r', 'g', 'b', 'm'};
for i = c'
    f = H.cells.faces(H.cells.facePos(i):H.cells.facePos(i+1)-1, :);
    arrayfun(@(j)plotFaces(H, f(j), 'edgecolor', cols{j}, 'linewidth', 2), ...
        1:4);
end
title('Directions of Pebi-style radial grid faces')
legend('H', 'Face 1', 'Face 2', 'Face 3', 'Face 4')

%% Build the radial grid by 'tessellationGrid'
% The 'tessellationGrid' accepts the combination of points (geometrical 
% input) and connectivity list (topological input) and returns a valid
% grid structure. The connectivity list provides the indices of nodes
% specifying each polygonal cell, i.e. the 'tessellationGrid' restores the 
% complete grid-object sequence of nodes->faces->cells from nodes->cells.
% The cell / node numbering follows the sequence of connectivity list /
% generating points.
% The face numbering follows the combining numbering of two nodes 
% specifying the face.
% e.g. for a connectivity list of {[1, 2, 4, 5], [2, 3, 4]}:
% Face indices:    1  2  3  4  5  6
% Node1 indices:   1  1  2  2  2  3
% Node2 indices:   4  5  3  4  5  4
%
% Since the radial grid is structured, we can pick the connectivity list 
% from the Cartesian node distribution matrix. The nodes corresponding to 
% cell (i,j) is: {L(i,j), L(i+1,j), L(i+1,j+1), L(i,j+1)}
%
% Arrange the node to a Cartesian format   
%   ---->  Radial
%  |
%  V  Angular
%
%      1  2  3                      nR+1
%  1   *  *  *  *  *  *  ....  *  *  *
%  2   *  *  *  *  *  *  ....  *  *  *
%  3   *  *  *  *  *  *  ....  *  *  *
%      ......                   ......
%      *  *  *  *  *  *  ....  *  *  *
%  nA  *  *  *  *  *  *  ....  *  *  *
%  1   *  *  *  *  *  *  ....  *  *  *

% Make the connectivity list
np = size(p,1);
nd = reshape((1:np)', nA, nR+1);
nd = [nd; nd(1,:)];
i = repmat((1:nA)',  1, nR);
j = repmat((1:nR),  nA,  1);
t = arrayfun(@(i, j)[nd(i,j), nd(i+1,j), nd(i+1,j+1), nd(i,j+1)], i(:), ...
    j(:), 'UniformOutput', false);

% Build the radial grid
G  = tessellationGrid(p, t);

% Plot the grid and generating points
figure, hold on; axis equal off
plotGrid(G)
demoPlotLine(p, 'ks', 'r', 3)
title('Radial grid and generating points')

% The logical indices
[a, r] = ind2sub([nA, nR], (1:G.cells.num)');
figure
subplot(1,2,1), axis equal tight off
plotCellData(G, a)
title('Angular indices')
subplot(1,2,2), axis equal tight off
plotCellData(G, r)
title('Radial indices')

% According to the picking method of connectivity list, for every cell:
% Face 1:  Radial  -
% Face 2:  Angular + 
% Face 3:  Radial  +
% Face 4:  Angular - 
% If the points are generated from R+ to R-, the directions of face 1 and 3 
% will be R+ and R-
% If the points are generated in clockwise direction, the directions of 
% face 2 and 4 will be A- and A+
c = randperm(G.cells.num, 4);
figure, hold on; axis equal off
plotGrid(G, 'facecolor', 'none')
cols = {'r', 'g', 'b', 'm'};
for i = c
    f = G.cells.faces(G.cells.facePos(i):G.cells.facePos(i+1)-1, :);
    arrayfun(@(j)plotFaces(G, f(j), 'edgecolor', cols{j}, 'linewidth', 2), ...
        1:4);
end
title('Directions of radial grid faces')
legend('G', 'Face 1 (R-)', 'Face 2 (A+)', 'Face 3 (R+)', 'Face 4 (A-)')

%%
% The function 'buildRadialGrid' integrates above generation process, we 
% can call it to build a radial grid directly.
G1 = buildRadialGrid(p, nA, nR)
