% 2D Three-Phase Surfactant-Polymer Injection Case
clc
clear

mrstModule add ad-core ad-blackoil ad-eor ad-props ...
               deckformat mrst-gui

%% Set up model and initial conditions
% The data required for the example
% The following are all the files needed for this tutorial
% Two files are the data for the simulation of surfactant polymer flooding.
current_dir = fullfile(mrstPath('ad-eor'), 'examples', 'surfactant-polymer');
fn = fullfile(current_dir, 'SURFACTANTPOLYMER2D.DATA');
gravity reset on;

deck = readEclipseDeck(fn);
% The deck is using metric system, MRST uses SI unit internally
deck = convertDeckUnits(deck);

% Construct physical model, initial state and dynamic well controls.
[state0, model, schedule] = initEclipseProblemAD(deck);

% Add initial surfactant & polymer concentration
state0.cp   = zeros([model.G.cells.num, 1]);
state0.cs   = zeros([model.G.cells.num, 1]);
state0.cpmax =  zeros([model.G.cells.num, 1]);

%% Select nonlinear and linear solvers

% Using physically normalized residuals for non-linear convergence
% calcuation.
% model.nonlinearTolerance = 1e-2;
model.useCNVConvergence = true;
% model.usingShear=false;
% model.toleranceCNV = 1.e-2;
% model.toleranceMB = 1.e-3;
% model.FacilityModel.toleranceWellRate = 1e-3;

% Setting up the non-linear solver.
nonlinearsolver = NonLinearSolver();
nonlinearsolver.useRelaxation = true;
nonlinearsolver.maxTimestepCuts = 10;
%% Full case - for comparison when surfactant is working
close all, clear statesSP, clear wellSolsSP;
scheduleSP = schedule;
% scheduleSP.control(1).W(1).cp=0;
% scheduleSP.control(2).W(1).cp=0;
[wellSolsSP, statesSP, reportsSP] = ...
    simulateScheduleAD(state0, model, scheduleSP, ...
                    'NonLinearSolver', nonlinearsolver);

%%
clear gmodel statesGP wellSolsGP;
gmodelsp = GenericSurfactantPolymerModel(model.G, model.rock, model.fluid, deck, 'disgas', model.disgas, 'vapoil', model.vapoil);
% gmodelp.surfactant = false;
gmodelsp.nonlinearTolerance = 1e-2;
% gmodelp.toleranceCNV = 1.e-3;
% gmodelsp.usingShear=false;

gmodelsp = gmodelsp.validateModel();
gmodelp.FacilityModel.toleranceWellRate = 1e-2;
% schedule.control(2).W(1).cp = 0;
% 
scheduleGSP = schedule;
% gmodelsp.polymer=false;
% scheduleGSP.control(1).W(1).cp=0;
% scheduleGSP.control(2).W(1).cp=0;

[wellSolsGSP, statesGSP, reportsGSP] = ...
    simulateScheduleAD(state0, gmodelsp, scheduleGSP, ...
                    'NonLinearSolver', nonlinearsolver);
                
%%
% clear gmodel statesGP;
% gmodelsp = GenericSurfactantPolymerModel(model.G, model.rock, model.fluid, 'disgas', model.disgas, 'vapoil', model.vapoil);
% % gmodelsp.surfactant = false;
% gmodelsp.nonlinearTolerance = 1e-2;
% gmodelsp.usingShear=false;
% 
% gmodelsp = gmodelsp.validateModel();
% gmodelsp.FacilityModel.toleranceWellRate = 1e-3;
% % schedule.control(2).W(1).cp = 0;
% % 
% scheduleGSP = schedule;
% [wellSolsGSP, statesGSP, reportsGSP] = ...
%     simulateScheduleAD(state0, gmodelsp, schedule, ...
%                     'NonLinearSolver', nonlinearsolver);

%%
% props = gmodelp.validateModel();
% groups = props.getStateFunctionGroupings();
% for i = 1:numel(groups)
%     figure;
%     [h,g] = plotStateFunctionGroupings(groups{i},  'label', 'label');
% %    printStateFunctionGroupingTikz(g);
%     title(class(groups{i}));
% end

%%
clear bomodel;
bomodel = GenericBlackOilModel(model.G, model.rock, model.fluid, 'disgas', model.disgas, 'vapoil', model.vapoil);
bomodel.nonlinearTolerance = 1e-2;
[wellSolsBO, statesBO, reportsBO] = ...
    simulateScheduleAD(state0, bomodel, schedule, ...
                    'NonLinearSolver', nonlinearsolver);
%%
% props = bomodel.validateModel();
% groups = props.getStateFunctionGroupings();
% for i = 1:numel(groups)
%     figure;
%     plotStateFunctionGroupings(groups{i},  'label', 'label')
%     title(class(groups{i}));
% end

%% 
plotWellSols({wellSolsGSP, wellSolsSP, wellSolsBO}, 'datasetnames', {'GenericSP', 'OldModel', 'Blackoil'});