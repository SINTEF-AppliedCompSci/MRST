%% Ensemble of GPSNET models 
% GPSNET models constitutes of a network model representing the reservoir
% and well connections. The model consists of a set of 1D reservoir models
% between the wells. In the graph terminology, the wells are represented as
% edges and the 1D reservoirs connecting the wells are edges.
% 
% In this example we set up an ensemble of such models, and we run ensemble
% simulations with uncertain rock properties (each connection has
% homogeneous rock properties) and well production indices.

mrstModule add ad-core ad-blackoil mrst-gui ad-props ...
    example-suite incomp ensemble dd-models diagnostics

mrstVerbose off

%% Set the name of the base problem and ensemble size
% The base problem contains all properties that are common throughout the
% ensemble

baseProblemName = 'ensemble_base_problem_simple_gpsnet_model';
baseProblemOptions = {};

ensembleSize = 80;

%% Define where to store truth and ensemble simulations
directoryTruth = fullfile(mrstOutputDirectory(), ...
                          'historyMatching', 'truth', ...
                          baseProblemName);
                      
topDirectory = fullfile(mrstOutputDirectory(), ...
                        'historyMatching', 'tutorial', baseProblemName);

%% Create example
baseExample = MRSTExample(baseProblemName, ...
                          'deleteOldResults', false, ...
                          'plotNetwork', false);

                      
%% Create the full model so that we can generate observations
trueExample = MRSTExample(baseExample.options.fullExampleName, ...
                          'deleteOldResults', false);
trueProblem = trueExample.getPackedSimulationProblem('Directory', directoryTruth);

plotExample = false;
rerunTrueProblemFromScratch = false;
overwriteObservation = true;

if rerunTrueProblemFromScratch
    clearPackedSimulatorOutput(trueProblem);
end
simulatePackedProblem(trueProblem);
if plotExample
    [wellSols, states, reports] = getPackedSimulatorOutput(trueProblem);
    trueExample.plot(states);
end

% Read and store the observations in a QoI object
trueQoI = WellQoI('wellNames', {'P1', 'P2'}, ...
                  'fldname', {'qOs', 'qWs'});
trueQoI = trueQoI.validateQoI(trueProblem);

% Create a separate ResultHandler for the observations 
observationResultHandler = trueQoI.ResultHandler;
observationResultHandler.dataPrefix = 'observedQoI';

% Create a separate ResultHandler for the observations.
% Need to build a new ResultHandler from scratch, so that we do not
% overwrite the dataPrefix property of observationResultHandler
truthResultHandler = ResultHandler('dataPrefix', 'trueQoI', ...
                                   'writeToDisk', observationResultHandler.writeToDisk,...
                                   'dataDirectory', observationResultHandler.dataDirectory, ...
                                   'dataFolder', observationResultHandler.dataFolder, ...
                                   'cleardir', false);

% Define the observation uncertainty and perturb the observations
% accordingly
obsStdDev = 0.00004;
trueObservations  = trueQoI.getQoI(trueProblem); 
if numel(observationResultHandler.getValidIds) < 1 || overwriteObservation
    for w = 1:numel(trueQoI.wellNames)
        for f = 1:numel(trueQoI.fldname)
            perturbedObservations{w}{f} = trueObservations{w}{f} + randn(size(trueObservations{w}{f}))*obsStdDev;
        end
    end
    observationResultHandler{1} = {perturbedObservations};
end
if numel(truthResultHandler.getValidIds) < 1 || overwriteObservation
    truthResultHandler{1} = {trueObservations};
end


                      
%% Define samples that give different rock properties for each connection

rockData = cell(ensembleSize, 1);
for i = 1:ensembleSize
    tmpPoro = zeros(1,4);
    for j = 1:4
        while tmpPoro(j) > 0.4 || tmpPoro(j) < 0.2
            tmpPoro(j) = 0.3 + randn(1)*0.1;
        end
    end
    rockData{i}.poro = tmpPoro;
    rockData{i}.perm = rockData{i}.poro.^3.*(1e-5)^2./(0.81*72*(1-rockData{i}.poro).^2);
end

%% Create sample object
rockSamples = NetworkRockSamples('data', rockData, ...
                                 'connectionIndices', baseExample.options.connectionIndices)

%% Create QoI, with the observations
qoi = WellQoI(...
    'wellNames', {'P1', 'P2'}, ...
    'fldname', {'qOs', 'qWs'}, ...
    'observationResultHandler', observationResultHandler, ...
    'truthResultHandler', truthResultHandler, ...
    'observationCov', obsStdDev^2);


%% Create the ensemble
rockEnsemble = MRSTHistoryMatchingEnsemble(baseExample, rockSamples, qoi, ...
    'directory', fullfile(topDirectory, 'rock'), ...
    'simulationStrategy', 'parallel', ...
    'maxWorkers', 8, ...
    'reset', true, ...
    'verbose', true);

%% Run ensemble
rockEnsemble.simulateEnsembleMembers();

%% Get simulated observations
disp('simulated observations: ')
size(rockEnsemble.getEnsembleQoI())

%% Get the matrix of ensemble samples 
disp('Matrix of ensemble samples (parameters):')
size(rockEnsemble.getEnsembleSamples())


%% Do history matching
rockEnsemble.doHistoryMatching()

%% Run new ensemble
rockEnsemble.simulateEnsembleMembers();

%% Plot original and updated ensemble results
rockEnsemble.plotQoI('subplots', true, 'clearFigure', false, ...
    'cmapName', 'lines', ...
    'plotTruth', true, ...
    'legend', {'observations', 'truth', 'posterior mean', 'prior mean'});

                








%% History match with both rock properties and well indices


wellSampleData = cell(ensembleSize, 1);
for i = 1:ensembleSize
    wellSampleData{i}.WI = rand(1,4)*1e-11;
end
wellSamples = WellSamples('data', wellSampleData);

compSamples = CompositeSamples({rockSamples, wellSamples}, 'tensorProduct', false);


%% Create combo ensemble
compEnsemble = MRSTHistoryMatchingEnsemble(baseExample, compSamples, qoi, ...
    'directory', fullfile(topDirectory, 'compRockWI'), ...
    'simulationStrategy', 'parallel', ...
    'maxWorkers', 8, ...
    'reset', true, ...
    'verbose', true)

%% Run ensemble
compEnsemble.simulateEnsembleMembers();

%% Get simulated observations
disp('simulated observations: ')
size(compEnsemble.getEnsembleQoI())

%% Get the matrix of ensemble samples 
disp('Matrix of ensemble samples (parameters):')
size(compEnsemble.getEnsembleSamples())


%% Do history matching
compEnsemble.doHistoryMatching()

%% Run new ensemble
compEnsemble.simulateEnsembleMembers();


%% Plot original and updated ensemble results
compEnsemble.plotQoI('subplots', true, 'clearFigure', false, ...
    'cmapName', 'lines', ...
    'plotTruth', true, ...
    'legend', {'observations', 'truth', 'posterior mean', 'prior mean'});
