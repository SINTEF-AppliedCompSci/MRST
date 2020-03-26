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

%% Select nonlinear and linear solvers

% Using physically normalized residuals for non-linear convergence
% calcuation.
model.useCNVConvergence = true;

% Setting up the non-linear solver.
nonlinearsolver = NonLinearSolver();
nonlinearsolver.useRelaxation = true;
%% Full case - for comparison when surfactant is working
% close all
% [wellSolsSP, statesSP, reportsSP] = ...
%     simulateScheduleAD(state0, model, schedule, ...
%                     'NonLinearSolver', nonlinearsolver);
%%
scheduleP = schedule;
scheduleP.control(1).W(1).cs = 0;
% scheduleP.control(2).W(1).cp = 0;
model.usingShear=false;

[wellSolsP, statesP, reportsP] = ...
    simulateScheduleAD(state0, model, scheduleP, ...
                    'NonLinearSolver', nonlinearsolver);
%%
clear gmodel statesGP;
gmodel = GenericSurfactantPolymerModel(model.G, model.rock, model.fluid, 'disgas', model.disgas, 'vapoil', model.vapoil);
gmodel.surfactant = false;
gmodel.nonlinearTolerance = 1e-2;
gmodel.usingShear=false;

gmodel = gmodel.validateModel();
gmodel.FacilityModel.toleranceWellRate = 1e-3;
% schedule.control(2).W(1).cp = 0;

scheduleGP = schedule;
[wellSolsGP, statesGP, reportsGP] = ...
    simulateScheduleAD(state0, gmodel, schedule, ...
                    'NonLinearSolver', nonlinearsolver);

%%
% props = gmodel.validateModel();
% groups = props.getStateFunctionGroupings();
% for i = 4:numel(groups)
%     figure;
%     [h,g] = plotStateFunctionGroupings(groups{i},  'label', 'label');
% %    printStateFunctionGroupingTikz(g);
%     title(class(groups{i}));
% end

%%
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
plotWellSols({wellSolsGP, wellSolsP, wellSolsBO}, 'datasetnames', {'New', 'Old', 'NoPolymer'})




%%
% gmodel2 = GenericSurfactantPolymerModel(model.G, model.rock, model.fluid, 'disgas', model.disgas, 'vapoil', model.vapoil);
% gmodel2.surfactant = false;
% gmodel2.nonlinearTolerance = 1e-2;
% gmodel2.usingShear=false;
% % schedule.control(2).W(1).cp = 0;
%
% [wellSolsGP2, statesGP2, reportsGP2] = ...
%     simulateScheduleAD(state0, gmodel2, schedule, ...
%                     'NonLinearSolver', nonlinearsolver);
% %%
% plotWellSols({wellSolsGP2, wellSolsP, wellSolsBO}, 'datasetnames', {'New', 'Old', 'NoPolymer'})
