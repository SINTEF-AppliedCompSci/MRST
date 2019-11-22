%% Demonstration of build an individual radial grid
%
% The radial grid offers superior performance in describing the radial flow
% in the well vicinity. The mrst grid factory has provided the grid
% constructor 'pebi' which can build a radial grid from a set of radially
% distributed points (see 'showPEBI' in 'modules\book\examples\grids'). 
% However, this method has two restrictions:
% The generating points will not be the grid nodes, instead they are
% loacted inside the grid cells, which goes against 
%
%
% are loacted at the grid cell center, not the grid nodes (see the
% definition of the Voronoi diagram)
% 
% However, the indices of grid objects, i.e. cells, faces, and nodes, are 
% not arranged logically. The logical indices are important for the radial 
% grid when assigning the rock properties, computing the radial 
% transmissibility, and finding the well cells, etc. This example 
% demonstrates how to build an radial grid using the constructor
% 'tessellationGrid'.
%

%% Generate radial point set
% nA - Number of cells in angular direction
% nR - Number of cells in radial direction
% rW - The minimum radius
% rM - The maximum radius

[nA, nR, rW, rM] = deal(40, 10, 2, 10);
th = linspace(0, 2*pi, nA+1); th = th(1:end-1);
r = logspace(log10(rW), log10(rM), nR+1);
[R, TH] = meshgrid(r, th);
[px, py] = pol2cart(TH(:), R(:));
% The points obey the logical numbering (angularly cycle fastest, 
% then radially)
p = [px(:), py(:)];

%% Build the radial grid by 'pebi'
% Note 'pebi' requires the continuous region, i.e. no 'hole' inside the
% region. So, the central points should be added into the generation point
% set.
pP = [[0, 0]; p];
H = pebi(triangleGrid(pP));
H = computeGeometry(H);

% Plot the grid and points. The points are the 
figure, hold on; axis equal tight off
plotGrid(H)
demoPlotLine(pP, 'ks', 'r', 2)

%%
min(H.faces.areas)

D = euclideanDistance(H.nodes.coords, H.nodes.coords);
D = sort(D(:)); D = D(D~=0);
min(D)
%%
figure, hold on; axis equal tight off
plotCellData(H, (1:H.cells.num)')
demoPlotLine(H.cells.centroids, 'ks-', 'r', 2)
%%
unique(diff(H.cells.facePos))
