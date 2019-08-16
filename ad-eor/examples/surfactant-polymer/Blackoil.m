
% clc
% clear

mrstModule add ad-core ad-blackoil ad-eor ad-props ...
               deckformat mrst-gui

current_dir = fileparts(mfilename('fullpath'));
fn = fullfile(current_dir, 'Blackoil.DATA');
% gravity reset on;
% 
% deck = readEclipseDeck(fn);
% deck = convertDeckUnits(deck);
% [state0, model, schedule] = initEclipseProblemAD(deck);
%%
gravity on

deck = readEclipseDeck(fn);
deck = convertDeckUnits(deck);

fluid = initDeckADIFluid(deck);
G = initEclipseGrid(deck);
G = computeGeometry(G);
rock  = initEclipseRock(deck);
rock  = compressRock(rock, G.cells.indexMap);
model = ThreePhaseBlackOilModel(G, rock, fluid, ...
                                                  'inputdata', deck, ...
                                                  'extraStateOutput', true);
                                              
schedule = convertDeckScheduleToMRST(model, deck);
state0B = initStateDeck(model,deck);
%% Select nonlinear and linear solvers

% Using physically normalized residuals for non-linear convergence
% calcuation.
model.useCNVConvergence = true;

% Setting up the non-linear solver.
% Since the size of the current prolbem is relatively small, the direct
% linear solver provided by MATLAB is efficient enough, so we are not using
% CPR-preconditioned linear sover by specifying 'useCPR' to be false.
% For bigger problem, you should try to activate CPR-preconditioned linear
% solver by specifying 'useCPR' to be true. By doing that, AGMG algebraic
% multigrid solver will be used if it is present. For much larger cases,
% the AGMG will improve the solution speed significantly and be required.
nonlinearsolver = getNonLinearSolver(model, 'DynamicTimesteps', false, ...
                                     'useCPR', false);
nonlinearsolver.useRelaxation = true;

%% Run the schedule with plotting function
% Once a system has been created it is trivial to run the schedule. Any
% options such as maximum non-linear iterations and tolerance can be set in
% the system struct.

% The AD-solvers allow for dyanmic plotting during the simulation process.
% We set up the following function to plot the evolution of the related
% variables (s:2 means oil saturation by default), the change of the well
% curves, and the a panel showing simulation progress. You can customize
% the function based on your own preference.
fn = getPlotAfterStep(state0B, model, schedule, ...
    'plotWell', true, 'plot1d', true,'axis tight', false, ...
    'field', 's:2');
[wellSolsB, statesB, reportB] = ...
    simulateScheduleAD(state0B, model, schedule, ...
                    'NonLinearSolver', nonlinearsolver, 'afterStepFn', fn);