%% CO2 injection with mobility control
% The co2-foam module provides simulation of CO2 injection into aquifers
% with mobility control using surfactant. The surfactant can dissolve in
% both brine and CO2 phases, according to a partitioning coefficient, and
% adsorbe to reservoir rock. Foam is assumed to be generated whereever the
% surfactant concentration is large enough. Foam modifies the CO2 mobility
% according to surfactant concentration, brine saturation and (optional)
% the gas velocity.

%% Add modules
mrstModule add co2-foam
mrstModule add ad-core ad-props ad-blackoil ad-eor
mrstModule add test-suite
mrstModule add mrst-gui

%% Set up test case
% We consider a box model of 1400 x 1400 x 100 m, discretized with 14 x 14
% x 7 grid blocks. We inject CO2 with foam-generating surfactant for the
% first four years, followed by 120 years of injection without surfactant.
cartDims = [14,14,7];
setupFoam1 = TestCase('qfs_foam', 'name', 'foam-no-ads', 'cartDims', cartDims);

%% Inspect model
% Surfactant properties are defined using the fluid object
model = setupFoam1.model.validateModel();
disp(model.fluid);
disp(model.fluid.foam);

% The effect of foam is modelled as a multiplier of the standard black-oil
% gas mobility, implemeted in the state function MultipliedMobility. The
% foam can partition into the gas phase, water phase, or both, and the
% partitioning is computed using the ConcentrationsPartitioning state
% function.
disp(model.FlowPropertyFunctions);

%% Simulate case
% Pack up the model problem and send to solver.
problemFoam1 = setupFoam1.getPackedSimulationProblem();
simulatePackedProblem(problemFoam1, 'restartStep', 1);

%% Inspect results
[wellSolsFoam1, statesFoam1] = getPackedSimulatorOutput(problemFoam1);
setupFoam1.plot(statesFoam1);

%% Simulate with adsorption
setupFoam2 = TestCase('qfs_foam', 'name', 'foam-ads', 'cartDims', cartDims, 'useAds', true);
problemFoam2 = setupFoam2.getPackedSimulationProblem();
simulatePackedProblem(problemFoam2, 'restartStep', 1);

%% Inspect results
[wellSolsFoam2, statesFoam2] = getPackedSimulatorOutput(problemFoam2);
setupFoam2.plot(statesFoam2);

%% Compare with no foam
% Simulate case with no surfactant injected (concentration set to zero).
setupNoFoam = TestCase('qfs_foam', 'name', 'no-foam', 'cartDims', cartDims);
setupNoFoam.schedule.control(1).W(1).cs = 0;
problemNoFoam = setupNoFoam.getPackedSimulationProblem();
simulatePackedProblem(problemNoFoam, 'restartStep', 1);

%% Inspect results
[wellSolsNoFoam, statesNoFoam] = getPackedSimulatorOutput(problemNoFoam);
setupNoFoam.plot(statesNoFoam);

%% Compare results
plotWellSols({wellSolsFoam1, wellSolsFoam2, wellSolsNoFoam}, ...
    setupFoam1.schedule.step.val, 'DatasetNames', ...
    {'No adsorption', 'Adsorption', 'No foam'});