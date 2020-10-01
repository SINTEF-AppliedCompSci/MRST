%% Ensemble of 3D reservoirs
% This is a minimal example (and a sanity test) for creating an ensemble of
% 3D reservoirs with two phases (oil and water).

mrstModule add ad-core ad-blackoil mrst-gui ad-props ...
    example-suite incomp

mrstVerbose off

%% Choose a base problem and ensemble size
% The base problem contains all properties that are common throughout the
% ensemble

baseProblemName = 'ensemble_base_problem_3d_reservoir';

ensembleSize = 20;


%% Run the base problem if we wish
% Simulate and plot it for illustration:
simulateAndPlotExample = true;
rerunBaseProblemFromScratch = false;
if simulateAndPlotExample
    baseExample = MRSTExample(baseProblemName);
    problem = baseExample.getPackedSimulationProblem();
    if rerunBaseProblemFromScratch
        clearPackedSimulatorOutput(problem);
    end
    simulatePackedProblem(problem);

    [wellSols, states, reports] = getPackedSimulatorOutput(problem);
    baseExample.plot(states);
end

%% Select and populate samples for stochastic configurations class

numCells = [10 6 4];

configData = cell(ensembleSize, 1);
for i = 1:ensembleSize
    configData{i}.porosity = gaussianField(numCells, [0.2 0.4]); 
    configData{i}.permeability = configData{i}.porosity.^3.*(1e-5)^2./(0.81*72*(1-configData{i}.porosity).^2);
end

stoch = StochasticRockConfigurations(configData);

%% Select quantity of interest class

qoi = WellProductionQuantityOfInterest(...
    'wellNames', {'P1', 'P2'}, ...
    'cumulative', false, ...
    'numTimesteps', []);


%% Create the ensemble

ensemble = MRSTEnsemble(baseProblemName, stoch, qoi, baseProblemOptions{:});

