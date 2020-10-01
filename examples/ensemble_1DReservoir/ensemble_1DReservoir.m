%% Ensemble of 1D reservoirs
% This is a minimal example for creating an ensemble of 1D reservoirs with
% two phases (oil and water).

mrstModule add ad-core ad-blackoil mrst-gui ad-props ...
    example-suite incomp

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
simulateAndPlotExample = true;
rerunBaseProblemFromScratch = true;
if simulateAndPlotExample
    baseExample = MRSTExample(baseProblemName, baseProblemOptions{:});
    problem = baseExample.getPackedSimulationProblem();
    if rerunBaseProblemFromScratch
        clearPackedSimulatorOutput(problem);
    end
    simulatePackedProblem(problem);

    [wellSols, states, reports] = getPackedSimulatorOutput(problem);
    baseExample.plot(states);
end


%% Select and populate samples for stochastic configurations class

configData = cell(ensembleSize, 1);
for i = 1:ensembleSize
    configData{i}.porosity = gaussianField([numCells 1 1], [0.2 0.4]); 
    configData{i}.permeability = configData{i}.porosity.^3.*(1e-5)^2./(0.81*72*(1-configData{i}.porosity).^2);
end

stoch = StochasticRockConfigurations(configData);

%% Select quantity of interest class

qoiName = 'WellProductionQuantityOfInterest';
qoiOptions = {'cumulative', true, 'numTimesteps', 10};


