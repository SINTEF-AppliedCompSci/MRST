%% 2D Three-Phase Surfactant-Polymer Injection Case
% Test of the generic and legacy solvers for EOR
mrstModule add ad-core ad-blackoil ad-eor ad-props ...
               deckformat mrst-gui

%% Set up model and initial conditions
% The data required for the example
% The following are all the files needed for this tutorial
% Two files are the data for the simulation of surfactant polymer flooding.
fn = fullfile(mrstPath('ad-eor'), 'examples', 'surfactant-polymer', 'SURFACTANTPOLYMER2D.DATA');
% Construct physical model, initial state and dynamic well controls.
[state0, model, schedule, nls] = initEclipseProblemAD(fn, 'timestepstrategy', 'none');
arg = {model.G, model.rock, model.fluid, ...
      'disgas', model.disgas, 'vapoil', model.vapoil,...
      'inputdata', model.inputdata};
%% 
lmodelsp = ThreePhaseSurfactantPolymerModel(arg{:});
close all, clear statesSP, clear wellSolsSP;
scheduleSP = schedule;
[wellSolsSP, statesSP, reportsSP] = ...
    simulateScheduleAD(state0, lmodelsp, scheduleSP, ...
                    'NonLinearSolver', nls);

%%
clear gmodel statesGP wellSolsGP;
gmodelsp = GenericSurfactantPolymerModel(arg{:});
scheduleGSP = schedule;

[wellSolsGSP, statesGSP, reportsGSP] = ...
    simulateScheduleAD(state0, gmodelsp, scheduleGSP, ...
                    'NonLinearSolver', nls);

%%
clear bomodel;
bomodel = GenericBlackOilModel(arg{:});
[wellSolsBO, statesBO, reportsBO] = ...
    simulateScheduleAD(state0, bomodel, schedule, ...
                    'NonLinearSolver', nls);

%% Compare results
plotWellSols({wellSolsGSP, wellSolsSP, wellSolsBO}, 'datasetnames', {'GenericSP', 'OldModel', 'Blackoil'});
%% Plot BO model functions
plotStateFunctionGroupings(bomodel)
title('Black-oil')
plotStateFunctionGroupings(gmodelsp)
title('EOR')
