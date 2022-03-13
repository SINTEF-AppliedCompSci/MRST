%% Effects of scales structures on trapping capacity
% The CO2 Storage Atlas grids are very coarse, even at the finest
% resolution supplied. As they cover vast scales and were meant for
% mapping, this is to be expected.
%
% To give an indication of the subscale morphology and potential traps that
% are not resolved in the atlas models, we will use a data model of
% Sleipner, which is the world's first subsea CO2 storage site and is
% located in the larger Utsira formation. The model was originally
% developed by Statoil has much higher spatial resolution. The model was
% publicly released by IEA and could be downloaded from their webpage by
% registred users:
% http://www.ieaghg.org/index.php?/20110329248/sleipner-benchmark-model.html
%
% Altogether, the example shows the importance of using fine-scale
% resolution when simulating long-term migration because coarsening may
% have a large impact on the trapped CO2 volumes.

mrstModule add co2lab
%% Load data and create grids
fprintf('Constructing Sleipner model...');
datapath = getDatasetPath('sleipner');
sleipner_deck = readGRDECL(fullfile(datapath, 'M9X1.grdecl'));

% Do mapaxis explicitly to get coinciding coordinate systems
ma = [436914 6475050 436914 6469150 440114 6469150];
coord = reshape(sleipner_deck.COORD,3,[])';
coord(:,1:2) = mapAxes(coord(:,1:2), ma);
coord = coord';
sleipner_deck.COORD=coord(:);

% Create top-surface grids for Sleipner
G_sleipner  = processGRDECL(sleipner_deck);
G_sleipner  = computeGeometry(G_sleipner);
Gt_sleipner = topSurfaceGrid(G_sleipner);
fprintf('done\n');

% Create top-surface grids for Utsira (uncoarsened)
fprintf('Constructing Utsira model...');
grdecl_utsira = getAtlasGrid('Utsirafm', 'coarsening', 1);
G_utsira = processGRDECL(grdecl_utsira{1});
Gt_utsira = topSurfaceGrid(G_utsira);
fprintf('done\n');

%% Trap analysis
res_sleipner = trapAnalysis(Gt_sleipner, true);
res_utsira = trapAnalysis(Gt_utsira, true);

%%
% We create a bounding box approximately equal to the fine Sleipner grid
% and use it to plot the corresponding area of the Utsira formation. The
% Sleipner grid is shown along with all local traps. As can be seen from
% the figure, what is a smooth surface in the coarse Utsira grid has
% several fine scale structures in the Sleipner grid, leading to several
% traps and potential rivers.

xs = Gt_sleipner.nodes.coords(:,1);
ys = Gt_sleipner.nodes.coords(:,2);

x = Gt_utsira.cells.centroids(:,1);
y = Gt_utsira.cells.centroids(:,2);
region = min(xs) < x & x < max(xs) & min(ys) < y & y < max(ys);

% Plot the grid
clf;
plotGrid(Gt_utsira, region, 'facec', [1 .6 .6])
plotGrid(Gt_sleipner, res_sleipner.traps == 0, 'facec', 'none', 'edgea',.2)
plotCellData(Gt_sleipner, res_sleipner.traps, res_sleipner.traps ~= 0, 'edgec', 'none')
view(-40, 50)
axis tight off
title('Utsira subset and Sleipner grid')

%% Estimate the unresolved local oscillations per area
% We find the total trap volume for Sleipner and divide it by the total
% area of the Sleipner case to find a rough estimate of the trap volume per
% area from small-scale oscillations.
%
% Note that we are always using volume in the geometrical sense: To find
% the amount of CO2 stored, both a porosity and a reference density of CO2
% is required.

trapvol_sleipner = sum(volumesOfTraps(Gt_sleipner, res_sleipner, []));
area_sleipner = sum(Gt_sleipner.cells.volumes);

finescaletraps = trapvol_sleipner/area_sleipner;
fprintf(['\nBy using the Atlas grid, approximately %2.5g liters of trapping\n'...
         'volume is not resolved per m^2 of area\n'], 1000*finescaletraps);

%% Extrapolate this estimate to the whole Utsira formation
trapvol_utsira = sum(volumesOfTraps(Gt_utsira, res_utsira, []));

area_utsira = sum(Gt_utsira.cells.volumes);
lost_volume = area_utsira*finescaletraps;
fprintf(['Total approximate unresolved trap volume for Utsira: %2.5g liters\n'...
        '(%1.2f%% of estimated large-scale trapped volume)\n'],...
        1000*lost_volume, 100*lost_volume./trapvol_utsira);

%% Get another estimate by removing global trends
% This estimate is obviously quite large as there may be global traps
% included in the fine Utsira grid that are counted twice. As we are
% primarily interested in structural traps that are smaller than what can
% be resolved by the coarse grid, we can obtain a more conservative
% estimate by removing the overall trends in the Utsira grid from the
% Sleipner grid and recomputing the traps. This is done by creating an
% interpolant from the Utsira top-surface grid and sampling the interpolant
% in the corresponding fine coordinates to obtain new z values.
%
% The three grids are plotted: Note how the adjusted grid (in green) has
% less curvature as it intersects the original grid, while having less traps
% volume. The largest trap is significantly reduced once the trend has been
% removed.

zinterp = TriScatteredInterp(Gt_utsira.cells.centroids(:, 1), ...
                             Gt_utsira.cells.centroids(:, 2), ...
                             Gt_utsira.cells.z(:));

Gt_adjusted = Gt_sleipner;

% Adjust z values by subtracting the interpolated z value and adding the
% average value in the area.
adjust = @(z, xy) z - zinterp(xy(:,1), xy(:,2)) + mean(Gt_utsira.cells.z(region));

Gt_adjusted.cells.z = adjust(Gt_adjusted.cells.z, Gt_adjusted.cells.centroids);
Gt_adjusted.faces.z = adjust(Gt_adjusted.faces.z, Gt_adjusted.faces.centroids);
Gt_adjusted.nodes.z = adjust(Gt_adjusted.nodes.z, Gt_adjusted.nodes.coords);

% Recompute geometry to get correct centroids
Gt_adjusted = computeGeometryVE_2D(Gt_adjusted);

res_adjusted = trapAnalysis(Gt_adjusted, true);

clf;
% Plot adjusted grid with traps
plotGrid(Gt_adjusted, res_adjusted.traps == 0, 'facea', .5, 'facec', 'green')
plotCellData(Gt_adjusted, res_adjusted.traps, res_adjusted.traps ~= 0, 'edgec','none')

% Plot original grid with traps
plotCellData(Gt_sleipner, res_sleipner.traps, res_sleipner.traps~=0, 'edgec','none')
plotGrid(Gt_sleipner, res_sleipner.traps == 0, 'facea', .3, 'facec', 'blue')

% Plot the Utsira grid showing the general trend being removed
plotGrid(Gt_utsira, region, 'facec', 'red')

view(-50, 40);
axis tight off

title('Adjusted versus original Sleipner grid')

%% Find new trapping volume
trapvol_adjusted = sum(volumesOfTraps(Gt_adjusted, res_adjusted, []));

finescaletraps = trapvol_adjusted/area_sleipner;
fprintf(['\nBy using the Atlas grid, approximately %2.5g liters of trapping\n'...
         'volume is unresolved per m^2 of area\n'], 1000*finescaletraps);

%% Extrapolate this estimate to the whole Utsira formation
lost_volume_adjusted = area_utsira*finescaletraps;
fprintf(['Total approximate unresolved trap volume for Utsira ' ...
        '(with global trends removed):\n'...
        '%2.5g liters (%1.2f%% of estimated large scale trapped volume)\n'],...
         1000*lost_volume_adjusted, 100*lost_volume_adjusted./trapvol_utsira);

%% Compare trap of local variations to global variations
% Shown as a pie chart, it is obvious that fine-scale variations in tje
% top-surface geometry contribute a significant volume to structural
% trapping compared to the large-scale structural traps. This shows the
% importance of using fine-scale resolution when simulating long-term
% migration because coarsening may have a large impact on the trapped CO2
% volumes.
figure
pie([lost_volume_adjusted trapvol_utsira])
legend({'Local variations', 'Global variations'})

%%
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
