%% Estimate the loss of information: Utsira formation
% The CO2 Storage Atlas grids are very coarse, even at the finest
% resolution supplied. As they cover vast scales and were meant for
% mapping, this is to be expected.
%
% To estimate the oscillations and potential traps lost in the coarsening
% process, the Sleipner dataset can be used. The Sleipner dataset is found
% at
% <http://www.ieaghg.org/index.php?/20110329248/sleipner-benchmark-model.html
% the IEAGHG website> (registration required) and contains fine scale data
% for the first subsea CO2 storage site. The Sleipner site is located in
% the larger Utsira formation and by comparing the local features we can
% see the level of detail which has been lost by coarsening.

mrstModule('add','co2lab', 'mex', 'coarsegrid', 'matlab_bgl', 'libgeometry', 'opm_gridprocessing');

sleipner_deck = readGRDECL(fullfile(mrstPath('co2lab'), 'data', 'sleipner', 'M9X1.grdecl'));

% Do mapaxis explicitly to get coinciding coordinate systems
%%{
ma = [436914 6475050 436914 6469150 440114 6469150];
coord = reshape(sleipner_deck.COORD,3,[])';
coord(:,1:2) = mapAxes(coord(:,1:2), ma);
coord = coord';
sleipner_deck.COORD=coord(:);
%}
% Create top surface grids for Sleipner and Utsira (uncoarsened)
try
   G_sleipner  = processgrid(sleipner_deck);
   G_sleipner  = mcomputeGeometry(G_sleipner);
catch
   G_sleipner  = processGRDECL(sleipner_deck);
   G_sleipner  = computeGeometry(G_sleipner);
end

Gt_sleipner = topSurfaceGrid(G_sleipner);

%%
grdecl_utsira = getAtlasGrid('Utsirafm', 'coarsening', 1);
try
   G_utsira = mprocessGRDECL(grdecl_utsira{1});
catch
   G_utsira = processGRDECL(grdecl_utsira{1});
end
Gt_utsira = topSurfaceGrid(G_utsira);

% Create trap analysis for Sleipner
res_sleipner = trapAnalysis(Gt_sleipner, true);
res_utsira = trapAnalysis(Gt_utsira, true);

%%
clf
p = get(gcf,'Position'); set(gcf,'Position', [p(1:2) 850 420]);
gu = Gt_utsira;
gu.nodes.coords = gu.nodes.coords/1e3;
gu.cells.centroids = gu.cells.centroids/1e3;
plotCellData(gu,gu.cells.z, 'EdgeColor','none');
view(-90,90); axis tight, box on
set(gca, 'Color',get(gcf,'Color'), 'FontSize',14);
h=colorbar; set(h,'FontSize',14);set(h,'XTick',.5,'XTickLabel','[m]')
xlabel('[km]'); ylabel('[km]');

[x,y,z]=deal(nan(gu.cartDims));
x(gu.cells.indexMap) = gu.cells.centroids(:,1);
y(gu.cells.indexMap) = gu.cells.centroids(:,2);
z(gu.cells.indexMap) = gu.cells.z;
hold on
contour(x,y,z,30,'k');

gs = Gt_sleipner;
gs.nodes.coords = gs.nodes.coords/1e3;
gs.cells.centroids = gs.cells.centroids/1e3;
bf = boundaryFaces(gs);
plotFaces(gs,bf,'EdgeColor','r','LineWidth',3);
hold off

%%
% We create a bounding box approximately equal to the fine Sleipner grid
% and use it to plot the corresponding area of the Utsira formation. The
% Sleipner grid is shown along with all local traps. As can be seen from
% the figure, what is a smooth surface in the coarse Utsira grid has
% several fine scale structures in the Sleipner grid, leading to several
% traps and potential rivers.

xs = gs.nodes.coords(:,1);
ys = gs.nodes.coords(:,2);

x = gu.cells.centroids(:,1);
y = gu.cells.centroids(:,2);
region = min(xs) < x & x < max(xs) & min(ys) < y & y < max(ys);

% Plot the grid
clf; set(gcf,'Position',p);
plotGrid(gu, region, 'facec', [.75 .5 .5])
plotGrid(gs, res_sleipner.traps == 0, 'facec', 'none','EdgeAlpha',.3)
plotCellData(gs, res_sleipner.traps, res_sleipner.traps ~= 0,'EdgeAlpha',.5)
view(-40, 50)
axis tight, box on

%% Estimate the lost local oscillations per area
% We find the total trap volume for Sleipner and divide it by the total
% area of the Sleipner case to find a rough estimate of the trap volume per
% area from small scale oscillations.
%
% Note that we are always using volume in the geometrical sense: To find
% the amount of CO2 stored both a porosity and a reference density of CO2
% is required.

trapvol_sleipner = sum(volumesOfTraps(Gt_sleipner, res_sleipner, []));
area_sleipner = sum(Gt_sleipner.cells.volumes);

finescaletraps = trapvol_sleipner/area_sleipner;
fprintf(['\nBy using the Atlas grid, approximately %2.5g liters of trapping\n'...
         'volume is lost per m^2 of area\n'], 1000*finescaletraps);

%% Extrapolate this estimate to the whole Utsira formation
trapvol_utsira = sum(volumesOfTraps(Gt_utsira, res_utsira, []));

area_utsira = sum(Gt_utsira.cells.volumes);
lost_volume = area_utsira*finescaletraps;
fprintf(['\nTotal approximate lost trap volume for Utsira: %2.5g m^3\n'...
        '(%1.2f%% of estimated large scale trapped volume)\n'],...
          lost_volume, 100*lost_volume./trapvol_utsira);

%% Get another estimate by removing global trends
% This estimate is obviously quite large as there may be global traps
% included in the fine Utsira grid which are then counted twice. As we are
% primarily interested in structural traps which are smaller than the
% coarse grid scale, we can obtain a more conservative estimate by removing
% the overall trends from the Utsira grid from the Sleipner grid and
% recomputing the traps.
%
% This is done by creating an interpolant from the Utsira top surface grid
% and sampling the interpolant in the corresponding fine coordinates to
% obtain new z values.
%
% The three grids are plotted: Note how the adjusted grid (in green) has
% less curvature as it intersects the original grid while having less traps
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

%%
clf;
% Plot original grid with traps
subplot(1,2,1);
%plotGrid(Gt_sleipner, res_sleipner.traps == 0, 'facecolor', 'none')
plotCellData(gs, res_sleipner.traps, res_sleipner.traps~=0,'EdgeColor','none')
plotFaces(gs,bf,'EdgeColor','k','LineWidth',3);
axis tight off;
title('Original', 'FontSize', 14)

% Plot adjusted grid with traps
subplot(1,2,2);
%plotGrid(Gt_adjusted, res_adjusted.traps == 0, 'facecolor', 'none')
plotCellData(gs, res_adjusted.traps, res_adjusted.traps ~= 0,'EdgeColor','none')
plotFaces(gs,bf,'EdgeColor','k','LineWidth',3);
axis tight off
title('Adjusted', 'FontSize', 14);

%set(gcf, 'color', [1 1 1]);

%% Find new trapping volume
trapvol_adjusted = sum(volumesOfTraps(Gt_adjusted, res_adjusted, []));

finescaletraps = trapvol_adjusted/area_sleipner;
fprintf(['\nBy using the Atlas grid, approximately %2.5g liters of trapping\n'...
         'volume is lost per m^2 of area\n'], 1000*finescaletraps);

%% Extrapolate this estimate to the whole Utsira formation
lost_volume_adjusted = area_utsira*finescaletraps;
fprintf(['\nTotal approximate lost trap volume for Utsira ' ...
        '(with global trends removed): '...
        '\n%2.5g m^3 (%1.2f%% of estimated large scale trapped volume)\n'],...
         lost_volume_adjusted, 100*lost_volume_adjusted./trapvol_utsira);