%% Multiscale solver applied to high-resolution bed model
% Pinchouts will create unstructured non-neighboring connections and hence
% be one of the principal gridding challenges for real-life reservoir
% models. To exemplify, we consider a highly detailed, core-scale model of
% realistic bedding structures. Such models are used as input to derive
% directional permeability for a given lithofacies and identify net pay
% below the level of petrophysical log resolution. The bedding structure is
% consists of six different rock types and is realized on a 30x30x333
% corner-point grid. Almost all cells are affected in some degree by pinch:
% Although the model has around 300.000 cells initially, over 2/3 will be
% inactive due to significant erosion, giving a fine grid of approximate
% dimensions 30x30x100. The volumes of the cells, as well as the areas of
% the (vertical) faces, vary almost four orders of magnitudes. We partition
% the fine grid into 6x6x5 coarse blocks and apply Dirichlet boundary
% conditions imposing a pressure drop from west to east. 
%
% The model not only has a difficult geometry, but also has a large number
% of low-permeable shale layers pinched between the other high-permeable
% layers. Impermeable regions like this are known to pose monotonicity
% problems for the MsFV method. Here, we demonstrate that this is not the
% case for the MsRSB method.

mrstModule add mrst-gui incomp coarsegrid msrsb

gravity off
%% Load model
% This is one of the standard data sets that comes along with MRST and can
% be downloaded automatically from our webpages.
pth = getDatasetPath('bedmodel2');

grdecl = readGRDECL(fullfile(pth, 'BedModel2.grdecl'));
usys   = getUnitSystem('METRIC');
grdecl = convertInputUnits(grdecl, usys);

G = processGRDECL(grdecl);
G = computeGeometry(G);

rock = grdecl2Rock(grdecl, G.cells.indexMap);

%% Plot the petrophysical data
% We launch two separate windows so that you can compare the petrophysical
% data and the corresponding rock types
figure;
plotToolbar(G, rock);
axis tight off
view(30, 30);

figure;
plotToolbar(G, grdecl.SATNUM(G.cells.indexMap));
axis tight off
view(30, 30);

%% Set boundary conditions
% Pressure drop from left to right and no-flow conditions elsewhere
xf = G.faces.centroids(:, 1);
left = find(abs(xf - min(xf)) < 1e-4);
right = find(abs(xf - max(xf)) < 1e-4);

bc = [];
bc = addBC(bc, left, 'pressure', 1);
bc = addBC(bc, right, 'pressure', 0);


figure;
plotGrid(G, 'facecolor', 'none');
plotFaces(G, left, 'facecolor', 'r')
plotFaces(G, right, 'facecolor', 'b')
view(30, 30)

%% Single-phase fluid model
fluid = initSingleFluid('rho', 1, 'mu', 1);
T = computeTrans(G, rock);
state = initResSol(G, 0);

%% Fine-scale reference solution
ref = incompTPFA(state, G, T, fluid, 'bc', bc);

%% Make coarse grid structure
% We partition the grid into 6x6x5 semiregular blocks by using the cell
% centroids to sample a partition vector from a 6x6x5 Cartesian
% representation of the bounding box. We then construct support regions.
p = sampleFromBox(G,reshape(1:6*6*5,[6,6,5]));
% p = partitionMETIS(G, T, 200);
CG = generateCoarseGrid(G, p);
CG = coarsenGeometry(CG);
CG = storeInteractionRegion(CG, 'edgeBoundaryCenters', true, 'adjustCenters', true);

%% Compute basis functions
A = getIncomp1PhMatrix(G, T, state, fluid);
basis = getMultiscaleBasis(CG, A, 'type', 'msrsb', 'iterations', 150);

%% Compute multiscale solution
ms = incompMultiscale(state, CG, T, fluid, basis, 'bc', bc);

%% Compare fine-scale and multiscale solutions
figure;
plotToolbar(G, ms), view(3)
title('MS')

figure;
plotToolbar(G, ref), view(3)
title('TPFA')

%% Copyright notice

% <html>
% <p><font size="-1">
% Copyright 2009-2018 SINTEF ICT, Applied Mathematics.
% </font></p>
% <p><font size="-1">
% This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).
% </font></p>
% <p><font size="-1">
% MRST is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% </font></p>
% <p><font size="-1">
% MRST is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% </font></p>
% <p><font size="-1">
% You should have received a copy of the GNU General Public License
% along with MRST.  If not, see
% <a href="http://www.gnu.org/licenses/">http://www.gnu.org/licenses</a>.
% </font></p>
% </html>
