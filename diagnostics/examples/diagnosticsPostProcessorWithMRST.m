%% Diagnostics GUI using MRST simulation output
%
% This example shows how to use the diagnostics post processing GUI and
% visualise results from MRST simulation output.
% First we run the Egg model in MRST using the packSimulationProblem
% setup. Then we load the results into the GUI which calculates the
% diagnostics and displays them interactively.
%
% See also PostProcessDiagnosticsECLIPSE.m
%
% Apr. 2019

mrstModule add ad-blackoil ad-core deckformat mrst-gui ad-props diagnostics
deck = getDeckEGG('realization',1);
gravity reset on

%% Setup model, schedule and initial state
G = initEclipseGrid(deck);
G = extractSubgrid(G, logical(deck.GRID.ACTNUM));
[state0, model, schedule, nonlinear] = initEclipseProblemAD(deck, 'G', G);
%% Define/pack simulation problem
problem = packSimulationProblem(state0, model, schedule, 'egg_model_FlowDiagnostics','NonLinearSolver',nonlinear);


%% Simulate problem
% If simulation is aborted, a new call to simulatePackedProblem will continue 
% the simulation. If simulation is complete, a new call will do nothing:
simulatePackedProblem(problem);

%% states/wellSols are accessed through handles
[ws, states, reports] = getPackedSimulatorOutput(problem);

%% Run PostProcessDiagnostics
% PostProcessDiagnosticsMRST takes as input a packed problem, packed using
% the packSimulationProblem function.
% Diagnostics are calculated and displayed interactively.

PostProcessDiagnosticsMRST(problem);
