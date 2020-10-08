%% Ensemble of 1D reservoirs
% This is a minimal example for creating an ensemble of 1D reservoirs with
% two phases (oil and water).

mrstModule add ad-core ad-blackoil mrst-gui ad-props ...
    example-suite incomp ensemble

mrstVerbose off

%% Choose a base problem and ensemble size
% The base problem contains all properties that are common throughout the
% ensemble

baseProblemName = 'ensemble_base_problem_1d_reservoir';
numCells = 30;
baseProblemOptions = {'ncells', numCells};

ensembleSize = 20;

%% Run the base problem if we wish
% Simulate and plot it for illustration:
simulateExample = true;
plotExample = false;
rerunBaseProblemFromScratch = false;

baseExample = MRSTExample(baseProblemName, baseProblemOptions{:});

if simulateExample
    problem = baseExample.getPackedSimulationProblem();
    if rerunBaseProblemFromScratch
        clearPackedSimulatorOutput(problem);
    end
    simulatePackedProblem(problem);

    [wellSols, states, reports] = getPackedSimulatorOutput(problem);
    if plotExample
        baseExample.plot(states);
    end
end


%% Select and populate samples for stochastic configurations class

configData = cell(ensembleSize, 1);
for i = 1:ensembleSize
    configData{i}.poro = gaussianField(baseExample.model.G.cartDims, [0.2 0.4]); 
    configData{i}.perm = configData{i}.poro.^3.*(1e-5)^2./(0.81*72*(1-configData{i}.poro).^2);
end

samples = RockSamples('data', configData);

%% Select quantity of interest class

qoi = WellQoI('wellNames', {'P1'}, 'cumulative', false, 'numTimesteps', 10, ...
    'fldname', {'qOs', 'qWs'});


%% Create the ensemble

ensemble = MRSTEnsemble(baseExample, samples, qoi, ... 
    ... %'directory', uniqueDirectory, ...
    'simulationType', 'parallel', ...
    'maxWorkers', 8, ...
    'verbose', true, ...
    'deleteOldResults', true...
    );


%% Run ensemble
ensemble.simulateAllEnsembleMembers();

%% Plot results
ensemble.qoi.plotEnsemble(ensemble);

