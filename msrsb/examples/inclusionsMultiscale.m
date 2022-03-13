%% Multiscale solver with inactive cells / inclusions
% A tough test case for any multiscale solver is a grid with many inactive
% cells that make coarse grids highly challenging. This example
% demonstrates the robustness of MsRSB when faced with a large number of
% inactive cells.
%
% This example corresponds to Example 4.3.5 in
%  A multiscale restriction-smoothed basis method for high contrast porous
%  media represented on unstructured grids. J. Comput. Phys, Vol. 304, pp.
%  46-71, 2016. DOI: 10.1016/j.jcp.2015.10.010

mrstModule add  msrsb coarsegrid mrst-gui incomp

%%  Grid and petrophysics
% We have stored a image containing the inclusions
img = imread(fullfile(mrstPath('msrsb'), 'examples', 'inclusions.png'));
img = double(img(:,:,1));

dims = size(img);
pdims = [1 1]*kilo*meter;
G = cartGrid(dims, pdims);
active = find(img == 255);

% Extract active cells
G = extractSubgrid(G, active);
G = computeGeometry(G);

% Create a simple log-normal permeability field and map onto active subset
k = logNormLayers(G.cartDims);
perm = k(G.cells.indexMap)*darcy;
poro = 0.5;
rock = makeRock(G, perm, poro);
T = computeTrans(G, rock);

clf,
plotCellData(G,log10(rock.perm),'EdgeColor','none')

% Single-phase unit fluid model. Since we are solving a single-phase
% incompressible problem without gravity, the values here do not matter.
fluid = initSimpleFluid('rho', [1, 1], 'mu', [1, 1], 'n', [1 1]);

%% Solve fine-scale problem
% We inject one pore volume over ten years with wells in opposing corners.
time = 10*year;
Nt = 100;
pv = poreVolume(G, rock);

% Set up wells in each corner
injr = sum(pv)/(time);
W = [];
W = verticalWell(W, G, rock, 1, 1, [], 'type', 'rate', ...
    'val', injr, 'comp_i', [1, 0]);
W = verticalWell(W, G, rock, G.cartDims(1), G.cartDims(2), [],...
    'type', 'bhp' , 'val', 100*barsa, 'comp_i', [1, 0]);

% Compute reference
state0 = initResSol(G, 0, [0, 1]);
ref0 = incompTPFA(state0, G, T, fluid, 'MatrixOutput', true, 'Wells', W);

clf, plotCellData(G,ref0.pressure/barsa,'EdgeColor','none');

%% Coarse partition
simple = true;
if simple
    % Simple coarse grid
    cdims = [10, 10];
    coarsedims = ceil(G.cartDims./cdims);
    p = partitionUI(G, coarsedims);
else
    % Alternative: Use Metis (if installed)
    p = partitionMETIS(G, T, 200); %#ok<UNRCH>
end
p = processPartition(G, p);
p = compressPartition(p);

% Plot the grid, log10 of permeability and coarse grid
clf
plotCellData(G, log10(rock.perm), 'EdgeColor', 'none')
outlineCoarseGrid(G, p, 'k')
title('Grid, permeability and coarse grid')
axis equal tight
colorbar

%% Coarse grid and basis functions
% Generate partition vector
CG = generateCoarseGrid(G, p);
% Add centroids / geometry information on coarse grid
CG = coarsenGeometry(CG);
% Store the support of each cell (used by multiscale basis construction)
CG = storeInteractionRegionCart(CG);
A = getIncomp1PhMatrix(G, T);
% Compute basis and solve for multiscale flux/pressure.
basis = getMultiscaleBasis(CG, A, 'type', 'msrsb');
ms0 = incompMultiscale(state0, CG, T, fluid, basis, 'wells', W);

%% Solve the time loop, and plot saturations as they progress
dt = time/Nt;
tsolve = @(state) implicitTransport(state, G, dt, rock, fluid, 'wells', W);
ref = ref0;
ms = ms0;
figure(1);
df = get(0, 'DefaultFigurePosition');
set(gcf, 'Position', df.*[1, 1, 2, 1])
for j = 1:Nt
    fprintf('%d of %d \n', j, Nt);

    disp('Solving MsRSB')
    ms = [ms; tsolve(ms(end))];
    disp('Solving TPFA')
    ref = [ref; tsolve(ref(end))];
    
    figure(1); clf
    subplot(1, 2, 1)
    plotCellData(G, ref(end).s(:, 1), 'EdgeColor', 'none')
    axis equal tight off
    caxis([0, 1])
    title('Fine scale')
    
    subplot(1, 2, 2)
    plotCellData(G, ms(end).s(:, 1), 'EdgeColor', 'none')
    axis equal tight off
    caxis([0, 1])
    title('MsRSB')
    drawnow
end
%% Plot L1 norm of the error in saturation
% We plot the saturation error with respect to the number of pore volumes
% injected.
L1 = @(x, y) norm(x - y, 1)/norm(y, 1);
err = arrayfun(@(x, y) L1(x.s(:, 1), y.s(:, 1)), ms, ref);

figure;
plot((1:numel(err))./numel(err), err, 'k--', 'linewidth', 2)
xlabel('PVI')
ylabel('L_1 error of saturation')

%% Launch interactive plotting
plotter = @(x) plotToolbar(G, x, 'field', 's:2');
close all
cmap = jet();
h1 = figure; plotter(ref);
axis equal tight; colormap(cmap)
title('Fine scale')

c = caxis();
h2 = figure; plotter(ms);
axis equal tight; colormap(cmap)
title('MsRSB')
caxis(c);

%% Copyright notice

% <html>
% <p><font size="-1">
% Copyright 2009-2022 SINTEF Digital, Mathematics & Cybernetics.
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

