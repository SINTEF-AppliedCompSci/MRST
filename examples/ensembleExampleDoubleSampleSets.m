%% Ensemble of 3D reservoirs
% Whereas this also serves as a minimal example, it demonstrates how to use
% different sample classes and combine them so that an ensemble has
% uncertain parameters related to different properties in the model
% (uncertain rock and uncertain well properties).
% To show this we create a small ensemble of a 3D reservoir modelling two
% phase flow (oil and water) between two injectors and two producers.

mrstModule add ad-core ad-blackoil mrst-gui ad-props ...
    example-suite incomp ensemble 

mrstVerbose off




%% Choose a base problem and ensemble size
% The base problem contains all properties that are common throughout the
% ensemble, and here we have implemented our own super simple example (see 
% example_template.m for the MRSTExample template).

baseProblemName = 'ensemble_base_problem_3d_reservoir';
baseProblemOptions = {};


% Change these flags to investigate the baseExample
simulateExample = false;
plotSimulation = false;
rerunBaseProblemFromScratch = false;

baseExample = MRSTExample(baseProblemName);

if simulateExample
    problem = baseExample.getPackedSimulationProblem();
    if rerunBaseProblemFromScratch
        clearPackedSimulatorOutput(problem);
    end
    simulatePackedProblem(problem);

    [wellSols, states, reports] = getPackedSimulatorOutput(problem);
    if plotSimulation
        baseExample.plot(states);
    end
end

%% an ensemble of stochastic rock realizations 
% RockSamples is a superclass of BaseSamples, and thereby implements all
% the sample functionality required to map new rock data to the base
% problem. Here, we precompute 20 rock realizations and use these as
% stochastic samples.

ensembleSize = 20;

rockData = cell(ensembleSize, 1);
for i = 1:ensembleSize
    rockData{i}.poro = gaussianField(baseExample.model.G.cartDims, [0.2 0.4]); 
    rockData{i}.perm = rockData{i}.poro.^3.*(1e-5)^2./(0.81*72*(1-rockData{i}.poro).^2);
end

rockSamples = RockSamples('data', rockData);

%% Select quantity of interest class
% We here demonstrates how the WellQoI class can be used to store
% production information for two fields for both producers.

qoi = WellQoI(...
    'wellNames', {'P1', 'P2'}, ...
    'fldname', {'qOs', 'qWs'}, ...
    'cumulative', false);


%% Create ensemble
% Combining the baseExample, the rock samples and the well QoI.
% We will run the ensemble in parallel across 8 workers, and we choose to
% delete any existing simulation results and run the ensemble from scratch

rockEnsemble = MRSTEnsemble(baseExample, rockSamples, qoi, ...
    'simulationStrategy', 'parallel', ...
    'maxWorkers', 8, ...
    'reset', true, ...
    'verbose', true);

%% Run ensemble
rockEnsemble.simulateEnsembleMembers();

%% Plot results
rockEnsemble.plotQoI();




%-------------------------------------------------------------------------%
%% Create another ensemble using stochastic well indices.
% The class WellSamples are used in the same way as RockSamples, and are
% also a superclass of BaseSamples. Instead of rock properties, however,
% its data property now holds well properties. In this example, we have
% stochastic well production indices (WI).

wellSampleData = cell(ensembleSize, 1);
for i = 1:ensembleSize
    wellSampleData{i}.WI = rand(1,4)*1e-11;
end

wellSamples = WellSamples('data', wellSampleData);

%% Define new ensemble
wellEnsemble = MRSTEnsemble(baseExample, wellSamples, qoi, ...
    'simulationStrategy', 'parallel', ...
    'maxWorkers', 8, ...
    'reset', true, ...
    'verbose', true);

%% Simulate and plot
wellEnsemble.simulateEnsembleMembers();

%% Plot result 
wellEnsemble.plotQoI();

%-------------------------------------------------------------------------%
%% Create another ensemble that combines both sampling strategies
% Finally, we want to run an ensemble that has stochastic rock properties
% AND stochastic well properties. To achieve this, we have implemented the
% class WellRockSamples that is a superclass of both WellSamples and
% RockSamples. 
% The data property in this class expects to have a cell array of structs
% with fields rock and well, both these being structs with rock properties
% and well properties, respectively.

comboData = cell(ensembleSize, 1);
for i = 1:ensembleSize
    comboData{i}.rock = rockData{i};
    comboData{i}.well = wellSampleData{i};
end

comboSamples = WellRockSamples('data', comboData);

%% Define new ensemble
comboEnsemble = MRSTEnsemble(baseExample, comboSamples, qoi, ...
    'simulationStrategy', 'parallel', ...
    'maxWorkers', 8, ...
    'reset', true, ...
    'verbose', true);

%% Simulate and plot
comboEnsemble.simulateEnsembleMembers();

%%
comboEnsemble.plotQoI();

