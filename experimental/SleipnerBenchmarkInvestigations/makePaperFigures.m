%% Figure generation for Sleipner Draft Paper.


%% Figure 1
% Observed CO2 outlines (plumes) superimposed onto a grid (Gt) to show mis-match
plumes = getLayer9CO2plumeOutlines();
[ ~, Gt, ~, ~ ] = makeSleipnerModelGrid('modelName','IEAGHGmodel', 'refineLevel',-4, 'plotsOn',false);
[plumes, topsurface, topfit, hCO2] = makeSurfaceDataAndPlots(plumes, Gt, 'plotsOn',true);


%% Figure 2
% CO2 entry rates into layer 9
[opt, var, model, schedule, initState] = studySleipnerBenchmarkFUN('ratecase','SPE', 'runSimulation',false);
[opt, var, model, schedule, initState] = studySleipnerBenchmarkFUN('ratecase','original', 'runSimulation',false);


%% Figure 3 - 4
% Fig3: GHGT model and IEAGHG model grids
% Fig4: Top surface comparison and re-centered difference
inspectSleipnerGridModels


%% Figure 5





%% Figure 6



%% Figure 7




%% Figure 8




%% Figure X
% Sensitivities:

% first run a simulation using smodel, and get the CO2 height data
[opt, var, smodel, schedule, initState, wellSols, states, sim_report] = studySleipnerBenchmarkFUN('refineLevel',-6, 'num_years',10, 'useSensModel',true);
states = addHeightData(states, smodel.G, smodel.fluid);

% then assess match between the 'simulated' states.h just obtained and the
% 'observed' states.h (i.e., from another simulation or using observed
% plume height data wrt a given grid)
plumes_base = getLayer9CO2plumeOutlines();
[plumes_base, ~, ~, ~] = makeSurfaceDataAndPlots(plumes_base, smodel.G);
dobj_dz = studySleipnerSensitivitiesFUN( initState, smodel, schedule, wellSols, states, 'plumes_base',plumes_base );






