%% Demonstration of gluing the radial grid to Cartesian grid
%
% This example demonstrates how to build the hybrid grid by gluing the 
% radial grid in the near-well region to the Cartesian grid elsewhere in
% the reservoir. The use of the constructor 'tessellationGrid' makes it
% easy to connect the two subgrids: obtain the points and connectivity
% lists of the subgrids and then assemble them to get the global ones. In 
% addition, the common points at the grid boundaries are merged and mapped. 

clear
mrstModule add nwm

%% Build the Cartesian grid
GC = cartGrid([20, 20], [200, 200]);
GC = computeGeometry(GC);

% Define the well region by logical indices
% 10 <= i <= 14
% 10 <= j <= 14
ij = gridLogicalIndices(GC);
idxI = ij{1} >= 10 & ij{1} <= 14 & ij{2} >= 10 & ij{2} <= 14;

% Cells inside the region
cI = find(idxI); 
% Cells outside the region
cO = find(~idxI);

% Plot the cells
figure, hold on; axis equal off
plotGrid(GC, cO, 'facecolor', 'y')
plotGrid(GC, cI, 'facecolor', 'g')
legend('Cells outside the well reigon', 'Cells inside the well reigon')

% Extract the sorted boundary nodes of the reigon (in counterclockwise)
bn = extractBdyNodesCells(GC, cI);

%% Generate radial point set
% nA - Number of cells in angular direction
% nR - Number of cells in radial direction
% rW - The minimum radius
% rM - The maximum radius

% The angular dimension and grid angles are determined by the boundary
% nodes to conform with the Cartesian grid
nA = numel(bn);

% Place the well at the region center
pCI = GC.cells.centroids(cI, :);
pW  = 0.5*[min(pCI(:,1)) + max(pCI(:,1)), min(pCI(:,2)) + max(pCI(:,2))];
pbn  = GC.nodes.coords(bn, :);
pbn0 = bsxfun(@minus, pbn, pW);

% Compute the angles
th = cart2pol(pbn0(:,1), pbn0(:,2));

% Define the logarithmic grid radii
[nR, rW, rM] = deal(10, 0.2, 16);
r = logspace(log10(rW), log10(rM), nR+1);

% Get the radial grid points
[R, TH] = meshgrid(r, th);
[px, py] = pol2cart(TH(:), R(:));
pR = bsxfun(@plus, [px(:), py(:)], pW);

% The boundary points are the outermost angular points. The total radial
% dimension is then nR+1
pR = [pR; pbn];

%% Build the radial grid
[GR, tR] = buildRadialGrid(pR, nA, nR+1);

% Plot the radial grid
figure, hold on; axis equal off
plotGrid(GC, cO, 'facecolor', 'y')
plotGrid(GR, 'facecolor', 'g')
legend('Cells outside the well reigon', 'Cells inside the well reigon')

%% Glue the grids
% The constructor 'tessellationGrid' is used to glue the two grids. First,
% we should get the points and connectivity lists of the two subgrids. The 
% connectivity list provides the indices of nodes specifying each polygonal 
% cell, i.e. the nodes of cells.

% Points and connectivity list of the radial grid
pR = GR.nodes.coords;
% tR is obtained from 'buildRadialGrid' to keep consistency of indices

% Points and connectivity list of the Cartesian grid
% Remove the cells inside the well region first
[GC_Rem, ~, ~, mapn] = removeCells(GC, cI);
pC = GC_Rem.nodes.coords;

% Get the indices of boundary nodes in GC_Rem
bn = arrayfun(@(n)find(mapn == n), bn);

% Merge the common nodes (boundary nodes). The boundary node indices in pC 
% are replaced by the ones in pR
nNo  = (1:size(pC,1))';
% The non-boundary nodes indices
idx  = ~ismember(nNo, bn);
nNo(idx) = (1:nnz(idx))' + size(pR,1);
pC = pC(idx,:);
% We already know the bounday nodes in pR: the last nA nodes
nNo(bn)  = size(pR,1)+1 - (nA:-1:1)';

% Map the connectivity list of GC_Rem
[cnC, pos] = gridCellNodes(GC_Rem, (1:GC_Rem.cells.num));
cnC = nNo(cnC);
tC  = arrayfun(@(c)cnC(pos(c):pos(c+1)-1), (1:GC_Rem.cells.num)', ...
    'UniformOutput', false);

% Assemble the points and connectivity lists
p = [pR; pC];
% Sort the tC to same direction with tR (in clockwise)
tC = sortPtsClockWise(p, tC);
t = [tR; tC];

% Build the hybrid grid
G = tessellationGrid(p, t);
G = computeGeometry(G);

% Plot the hybrid grid
figure, hold on; axis equal off
plotGrid(G)

%%
% The function 'radCartHybridGrid' integrates above generation process, we 
% can call it to build a hybrid grid directly.
G1 = radCartHybridGrid(GC, cI, rW, rM, nR, pW)