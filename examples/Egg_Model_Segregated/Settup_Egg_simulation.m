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
%realization = [2];
[G, rock, fluid, deck] = setupEGG();

%% Changing the well cells for producer form 1:7 to 1:3
for i = 9:12
    deck.SCHEDULE.control.COMPDAT{i,5} = [3];
end

[state, model, schedule, nonlinear] = initEclipseProblemAD(deck, 'G', G, 'TimestepStrategy', 'none');

load('Egg_Rock.mat','rock')

model.getPhaseNames()

%% Adding oil top and water bottom
for i = 1:model.G.cells.num
    if model.G.cells.centroids(i,3) < 4014
        state.s(i,:) = [0 1];
    else % model.G.cells.centroids(i,3) >= 4010
        state.s(i,:) = [1 0];
    end
end

% plotGrid(G, 'FaceColor', 'none', 'EdgeAlpha', 0.1)
% plotWell(G,schedule.control.W)


problem = packSimulationProblem(state, model, schedule, 'EGG_realization_0_segregated', 'NonLinearSolver', nonlinear);

%%
[ok, status] = simulatePackedProblem(problem);


%% Run simulation
 [wellSols, states, reports] = getPackedSimulatorOutput(problem);
    
  d = PostProcessDiagnosticsMRST(problem);
