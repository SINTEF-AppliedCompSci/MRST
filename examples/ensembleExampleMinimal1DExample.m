%% Ensemble of 1D reservoirs
% This is a minimal example for creating an ensemble of 1D reservoirs with
% two phases (oil and water) driven by an injector and a producer well
% located at each end of the reservoir.

mrstModule add ad-core ad-blackoil mrst-gui ad-props ...
    example-suite incomp ensemble

mrstVerbose off

%% Choose a base problem and ensemble size
% The base problem contains all properties that are common throughout the
% ensemble, and here we have implemented our own super simple example (see 
% example_template.m for the MRSTExample template).

baseProblemName = 'ensemble_base_problem_1d_reservoir';
numCells = 30;
baseProblemOptions = {'ncells', numCells};

baseExample = MRSTExample(baseProblemName, baseProblemOptions{:});

% Change these flags to investigate the baseExample
simulateExample = false;
plotExample = false;
rerunBaseProblemFromScratch = false;

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


%% Create samples that represent the stochastic parameters in our ensemble

ensembleSize = 20;

configData = cell(ensembleSize, 1);
for i = 1:ensembleSize
    configData{i}.poro = gaussianField(baseExample.model.G.cartDims, [0.2 0.4]); 
    configData{i}.perm = configData{i}.poro.^3.*(1e-5)^2./(0.81*72*(1-configData{i}.poro).^2);
end

samples = RockSamples('data', configData)

%% Select quantity of interest class

qoi = WellQoI('wellNames', {'P1'}, 'cumulative', true, ...
    'fldname', {'qOs', 'qWs'})


%% Create the ensemble

ensemble = MRSTEnsemble(baseExample, samples, qoi, ... 
    'simulationStrategy', 'parallel', ...
    'maxWorkers', 8, ...
    'verbose', true, ...
    'reset', true...
    )


%% Run ensemble
ensemble.simulateEnsembleMembers();

%% Plot results
ensemble.plotQoI();
