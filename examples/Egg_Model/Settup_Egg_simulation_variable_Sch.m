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
realization = 1;
[G, rock, fluid, deck] = setupEGG('realization', realization);
[state, model, schedule, nonlinear] = initEclipseProblemAD(deck, 'G', G, 'TimestepStrategy', 'none');

load('Egg_Rock.mat','rock')

model.rock=rock;
model.getPhaseNames()

schedule.step.control(1:48)=1;

schedule.control(2) =schedule.control(1);
schedule.control(2).W(2).val =(60.5/79.5)* schedule.control(1).W(2).val;
schedule.control(2).W(4).val =(120.5/79.5)* schedule.control(1).W(4).val;
schedule.control(2).W(6).val =(60.5/79.5)* schedule.control(1).W(6).val;
schedule.control(2).W(8).val =(60.5/79.5)* schedule.control(1).W(8).val;

schedule.step.control(49:72)=2;

schedule.control(3) =schedule.control(1);
schedule.control(3).W(1).val =(100.5/79.5)* schedule.control(1).W(1).val;
schedule.control(3).W(3).val =(100.5/79.5)* schedule.control(1).W(3).val;
schedule.control(3).W(4).val =(60.5/79.5)* schedule.control(1).W(4).val;
schedule.control(3).W(5).val =(100.5/79.5)* schedule.control(1).W(5).val;
schedule.control(3).W(7).val =(100.5/79.5)* schedule.control(1).W(7).val;

schedule.step.control(73:96)=3;

schedule.control(4) =schedule.control(1);
schedule.control(4).W(9).val =(375/395)* schedule.control(1).W(9).val;
schedule.control(4).W(11).val =(400/395)* schedule.control(1).W(11).val;
schedule.control(4).W(12).val =(375/395)* schedule.control(1).W(12).val;

schedule.step.control(97:end)=4;

problem = packSimulationProblem(state, model, schedule, 'DECK_egg_Flownet_wells_variable', 'NonLinearSolver', nonlinear);

%%
[ok, status] = simulatePackedProblem(problem);


%% Run simulation
 [wellSols, states, reports] = getPackedSimulatorOutput(problem);
    
  d = PostProcessDiagnosticsMRST(problem);
