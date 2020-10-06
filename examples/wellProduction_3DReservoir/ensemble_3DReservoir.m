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
    'cumulative', false, ...
    'numTimesteps', []);

%% Define a unique folder in case we want to redo the simulations
uniqueDirectory = fullfile(mrstOutputDirectory(), 'ensemble', 'unique', datestr(datetime('now'), 'yymmdd_hhMMss'));


%% Create the ensemble
ensemble = MRSTEnsemble(baseExample, rockSamples, qoi, ...
    ... %'directory', uniqueDirectory, ...
    'simulationType', 'parallel', ...
    'maxWorkers', 8);

%% Run ensemble
ensemble.simulateAllEnsembleMembers();

%% Gather result
t = cumsum(baseExample.schedule.step.val);
p1_production = cell(ensembleSize, 1);
p2_production = cell(ensembleSize, 2);
meanProduction = zeros(numel(t), 2);
for i = 1:ensembleSize
    results = ensemble.qoi.ResultHandler{i}{1};
    meanProduction = meanProduction + results;
    p1_production{i} = results(:,1);
    p2_production{i} = results(:,2);
end
meanProduction = meanProduction/ensembleSize();

%% Plot results

figure
subplot(2,1,1);
hold on
for i = 1:ensembleSize
    plot(t, p1_production{i},  'color', [1 1 1]*0.6);
end
plot(t, meanProduction(:,1), 'color', 'red');
title(strcat('Oil production from ', ensemble.qoi.wellNames(1)));

subplot(2,1,2);
hold on 
for i = 1:ensembleSize
    plot(t, p2_production{i},  'color', [1 1 1]*0.6);
end
plot(t, meanProduction(:,2), 'color', 'red');
title(strcat('Oil production from ', ensemble.qoi.wellNames(2)));