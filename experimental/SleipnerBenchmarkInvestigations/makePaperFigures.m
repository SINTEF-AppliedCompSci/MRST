%% Figure generation for Sleipner Draft Paper.


%% Figure 1
% Observed CO2 outlines (plumes) superimposed onto a grid (Gt) to show mis-match
plumes = getLayer9CO2plumeOutlines();
[ ~, Gt, ~, ~ ] = makeSleipnerModelGrid('modelName','IEAGHGmodel', 'refineLevel',-4, 'plotsOn',false);
[plumes, topsurface, topfit, hCO2] = makeSurfaceDataAndPlots(plumes, Gt, 'plotsOn',true);


%% Figure 2
% CO2 entry rates into layer 9



%% Figure 3 - 4
% Fig3: GHGT model and IEAGHG model grids
% Fig4: Top surface comparison and re-centered difference
inspectSleipnerGridModels


%% Figure 5





%% Figure 6



%% Figure 7




%% Figure 8




%% Figure 10






