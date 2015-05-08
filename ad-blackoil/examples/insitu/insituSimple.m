%% Example demonstrating in-situ plotting capabilities in MRST-AD
% The AD-solvers allow for dynamic plotting during the simulation process.
% This example demonstrates this capability using "plotWellSols" and a
% panel showing simulation progress.
%
% We first set up a simulation model of SPE1 in the standard manner.
mrstModule add ad-core ad-blackoil ad-props ad-fi mrst-gui deckformat

[G, rock, fluid, deck, state0] = setupSPE1();
model = selectModelFromDeck(G, rock, fluid, deck);
schedule = convertDeckScheduleToMRST(G, model, rock, deck);

%% Set up plotting function
% We set up a function handle that takes in the current simulator
% variables, which will be run after each succesful control step. This
% function can also pause the simulation or change the maximum timestep as
% a proof of concept.
close all
fn = getPlotAfterStep(state0, model, schedule, 'plotWell', true,...
                                               'plotReservoir', false);
disp(fn)
%% Run the simulation with plotting function
linsolve = CPRSolverAD();
[wellSols, states, report] = simulateScheduleAD(state0, model, schedule, ...
    'Verbose', true, 'afterStepFn', fn, 'linearsolver', linsolve);
