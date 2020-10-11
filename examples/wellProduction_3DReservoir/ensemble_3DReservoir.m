%% Ensemble of 3D reservoirs
% This is a minimal example (and a sanity test) for creating an ensemble of
% 3D reservoirs with two phases (oil and water).

mrstModule add ad-core ad-blackoil mrst-gui ad-props ...
    example-suite incomp ensemble

mrstVerbose off




%% Choose a base problem and ensemble size
% The base problem contains all properties that are common throughout the
% ensemble

baseProblemName = 'ensemble_base_problem_3d_reservoir';
baseProblemOptions = {};

ensembleSize = 20;


%% Run the base problem if we wish
% Simulate and plot it for illustration:
simulateExample = true;
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

%% Select and populate samples for stochastic configurations class

configData = cell(ensembleSize, 1);
for i = 1:ensembleSize
    configData{i}.poro = gaussianField(baseExample.model.G.cartDims, [0.2 0.4]); 
    configData{i}.perm = configData{i}.poro.^3.*(1e-5)^2./(0.81*72*(1-configData{i}.poro).^2);
end

rockSamples = RockSamples('data', configData);

%% Select quantity of interest class

qoi = WellQoI(...
    'wellNames', {'P1', 'P2'}, ...
    'fldname', {'qOs', 'qWs'}, ...
    'cumulative', false, ...
    'numTimesteps', []);

%% Define a unique folder in case we want to redo the simulations
uniqueDirectory = fullfile(mrstOutputDirectory(), 'ensemble', 'unique', datestr(datetime('now'), 'yymmdd_hhMMss'));


%% Create the ensemble
ensemble = MRSTEnsemble(baseExample, rockSamples, qoi, ...
    ... %'directory', uniqueDirectory, ...
    'simulationType', 'parallel', ...
    'maxWorkers', 8, ...
    'deleteOldResults', true, ...
    'verbose', true);

%% Run ensemble
ensemble.simulateAllEnsembleMembers();

%% Plot results
ensemble.qoi.plotEnsemble(ensemble);





%% Create another ensemble using stochastic well indices
% ---------------------
wellSampleData = cell(ensembleSize, 1);
for i = 1:ensembleSize
    wellSampleData{i}.WI = rand(1,4)*1e-11;
end

wellSamples = WellSamples('data', wellSampleData);

%% Define new ensemble
ensemble = MRSTEnsemble(baseExample, wellSamples, qoi, ...
    ... %'directory', uniqueDirectory, ...
    'simulationType', 'parallel', ...
    'maxWorkers', 8, ...
    'deleteOldResults', true, ...
    'verbose', true);

%% Simulate and plot
ensemble.simulateAllEnsembleMembers();

%%
ensemble.qoi.plotEnsemble(ensemble);

