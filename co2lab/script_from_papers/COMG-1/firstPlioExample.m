%% Residual trapping in the Pliocenesand formation
% The actual sand body of the Pliocenesand formation lies too shallow to be
% considered as a real candidate for CO2 storage, but the model can be used
% as a synthetic test case if we increase its burial depth to, e.g., a
% thousand meters. The top surface has almost no fine-scale structure and
% thus allows for a very low percentage (0.02%) of structural trapping
% compared to the overall volume of the whole sand body. To store CO2, one
% can therefore not rely on structural trapping and should instead consider
% residual trapping.

mrstModule add co2lab incomp ad-core coarsegrid;
gravity reset on;

grdecl = getAtlasGrid('Pliocenesand');
G      = processGRDECL(grdecl{1});
G      = computeGeometry(G(1));
Gt     = topSurfaceGrid(G);
ta     = trapAnalysis(Gt, false);

%%
% The formation describes a ridge that ends up in a relatively large plain;
figure; plotGrid(G,'EdgeAlpha',.05);
axis tight; view(120,40);

%%
% To demonstrate that there is almost no potential for structural trapping,
% we just load the 3D viewer
h = interactiveTrapping(Gt, 'method', 'node', 'light', true, 'spillregions', true);
view(100,50);

%%
% Next, we will perform a VE simulation. This could, of course, have been
% launched from inside the interactive viewer, but to make the example as
% reproducible as possible, we launch it manually from the outside.
%
% Remark: change 'T_migration' to 50 years and 'Nm'to 10 to reproduce middle
% part of figure 8.
petrodata.avgperm = 1.2*darcy;
petrodata.avgporo = 0.25;

% Determine correct cell index for well position, based on geographical
% coordinate (latitude and longitude)
coord = [464328, 6646937]; 
dist2 = sum(bsxfun(@minus, Gt.cells.centroids, coord).^2, 2);
[~, wellcell] = min(dist2);

close all
migrateInjection(Gt, ta, petrodata, wellcell, ...
                 'amount',      10, ... % Mt/year
                 'T_injection', 50*year,   ...
                 'T_migration', 0*year, ...
                 'topPressure', 100*barsa, ...
                 'Ni',          25,  ... % time steps during injection
                 'Nm',          0,  ... % time steps during migration
                 'plotPanel',   true);

%%
% Rearrange the plot
chld = get(gcf,'Children');
view(chld(7),105,50);
delete(chld(5));
set(chld(1),'Position',[.835 .65 .14 .18])
set(chld(2),'Position',[.6 .425 .375 .15]);
set(chld(3),'Position',[.2 .425 .375 .15]);

%%
% Rerun for a longer time
petrodata.avgperm = 1.2*darcy;
petrodata.avgporo = 0.25;
[~,report] = migrateInjection(Gt, ta, petrodata, wellcell, ...
                 'amount',      10, ... % Mt/year
                 'T_injection', 50*year,   ...
                 'T_migration', 1450*year, ...
                 'topPressure', 100*barsa, ...
                 'Ni',          25,   ... % time steps during injection
                 'Nm',          145,  ... % time steps during migration
                 'plotPanel',   true, ...
                 'plotHist',    true);
