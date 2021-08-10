%% Example demonstrating the two-phase oil-water Egg model
% This example sets up and runs the Egg model using the two-phase AD
% solvers. 
%
% For details on the EggModel and the corresponding ensamble, see
% Jansen, J. D., et al. "The egg model–a geological ensemble for reservoir
% simulation." Geoscience Data Journal 1.2 (2014): 192-195.

mrstModule add ad-core ad-blackoil deckformat diagnostics

% Realizations can be set to 0 for base cae, or a number between 1 and 100
% for different permeabilities.sa

[G, rock, fluid, deck] = setupEGG('realization', realization);
[state, model, schedule, nonlinear] = initEclipseProblemAD(deck, 'G', G, 'TimestepStrategy', 'none');


model.getPhaseNames()


problem = packSimulationProblem(state, model, schedule,['EGG_realization_',num2str(realization),'_paper'], 'NonLinearSolver', nonlinear);

%%
[ok, status] = simulatePackedProblem(problem);


%% Run simulation
 [wellSols, states, reports] = getPackedSimulatorOutput(problem);
    
%  d = PostProcessDiagnosticsMRST(problem);
