%% Ensemble of GPSNET models 
% GPSNET models constitutes of a network model representing the reservoir
% and well connections. The model consists of a set of 1D reservoir models
% between the wells. In the graph terminology, the wells are represented as
% edges and the 1D reservoirs connecting the wells are edges.
% 
% In this example we set up an ensemble of such models, and we run ensemble
% simulations with uncertain operator properties with each connection
% having homogeneous operator properties throughout the connection.

mrstModule add ad-core ad-blackoil mrst-gui ad-props ...
    example-suite incomp ensemble dd-models diagnostics

mrstVerbose off

%% Set the name of the base problem and ensemble size
% The base problem contains all properties that are common throughout the
% ensemble

baseProblemName = 'ensemble_base_problem_simple_gpsnet_model';
baseProblemOptions = {};

ensembleSize = 70;
numConnections = 4; % Not a parameter, but hard-coded in the baseProblem

%% Define where to store truth and ensemble simulations
directoryTruth = fullfile(mrstOutputDirectory(), ...
                          'historyMatching', 'truth','network_state0', ...
                          baseProblemName);
                      
topDirectory = fullfile(mrstOutputDirectory(), ...
                        'historyMatching', 'tutorial', 'network_state0', baseProblemName);

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
trueQoI = WellQoIHM('wellNames', {'P1', 'P2'}, ...
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


                      
%% Define samples that give different transmissibilise for each connection

initSoData = cell(ensembleSize, 1);
for i = 1:ensembleSize
    initSoTmp = zeros(numConnections);
    while((min(initSoTmp(:)) < 0.2) || (max(initSoTmp(:)) > 1))
        initSoTmp = randn(numConnections);
        initSoTmp = initSoTmp - mean(initSoTmp(:));
        initSoMin = min(initSoTmp(:));
        initSoMax = max(initSoTmp(:));
        initSoTmp = initSoTmp.*(0.45/(initSoMax - initSoMin)) + 0.7;
    end
    
    initSoData{i}.initSo = initSoTmp;
end


%% Create sample object
initSoSample = NetworkState0SamplesHM('data', initSoData, ...
                                      'connectionIndices', baseExample.options.connectionIndices)

%% Create QoI, with the observations
qoi = WellQoIHM(...
    'wellNames', {'P1', 'P2'}, ...
    'fldname', {'qOs', 'qWs'}, ...
    'observationResultHandler', observationResultHandler, ...
    'truthResultHandler', truthResultHandler, ...
    'observationCov', obsStdDev^2);


%% Create the ensemble
initSoEnsemble = MRSTHistoryMatchingEnsemble(baseExample, initSoSample, qoi, ...
    'directory', fullfile(topDirectory, 'state0Network'), ...
    'simulationStrategy', 'parallel', ...
    'maxWorkers', 8, ...
    'reset', true, ...
    'verbose', true);

%% Run ensemble
initSoEnsemble.simulateEnsembleMembers();


%% Get simulated observations
disp('simulated observations: ')
size(initSoEnsemble.getEnsembleQoI())

%% Get the matrix of ensemble samples 
disp('Matrix of ensemble samples (parameters):')
size(initSoEnsemble.getEnsembleSamples())


%% Do history matching
initSoEnsemble.doHistoryMatching()

%% Run new ensemble
initSoEnsemble.simulateEnsembleMembers();

%% Plot original and updated ensemble results
initSoEnsemble.plotQoI('subplots', true, 'clearFigure', false, ...
    'cmapName', 'lines', ...
    'plotTruth', true, ...
    'legend', {'observations', 'truth', 'posterior mean', 'prior mean'});

%%                
crash







%% History match with porevolumes
FDmean = baseExample.model.operators.pv(baseExample.options.connectionIndices.cells{1}(1));
pvData = cell(ensembleSize, 1);
for i = 1:ensembleSize
    pvData{i}.pv = FDmean + 2000*randn(1, numConnections);
end

pvSamples = NetworkOperatorSamplesHM('data', pvData, ...
                                     'connectionIndices', baseExample.options.connectionIndices, ...
                                     'pvScale', baseExample.model.G.cells.volumes(1)/100)



%% Create combo ensemble
pvEnsemble = MRSTHistoryMatchingEnsemble(baseExample, pvSamples, qoi, ...
    'directory', fullfile(topDirectory, 'pv'), ...
    'simulationStrategy', 'parallel', ...
    'maxWorkers', 8, ...
    'reset', true, ...
    'verbose', true)

%% Run ensemble
pvEnsemble.simulateEnsembleMembers();


%% Get simulated observations
disp('simulated observations: ')
size(pvEnsemble.getEnsembleQoI())

%% Get the matrix of ensemble samples 
disp('Matrix of ensemble samples (parameters):')
size(pvEnsemble.getEnsembleSamples())


%% Do history matching
pvEnsemble.doHistoryMatching()

%% Run new ensemble
pvEnsemble.simulateEnsembleMembers();


%% Plot original and updated ensemble results
pvEnsemble.plotQoI('subplots', true, 'clearFigure', false, ...
    'cmapName', 'lines', ...
    'plotTruth', true, ...
    'legend', {'observations', 'truth', 'posterior mean', 'prior mean'});










%% History match with both porevolumes and transmissibilities
operatorData = cell(ensembleSize, 1);
for i = 1:ensembleSize
    operatorData{i}.pv = pvData{i}.pv;
    operatorData{i}.T = transData{i}.T;
end

operatorSamples = NetworkOperatorSamplesHM('data', operatorData, ...
                                     'connectionIndices', baseExample.options.connectionIndices, ...
                                     'pvScale', baseExample.model.G.cells.volumes(1)/100)


                                 
                                
%% Create combo ensemble
operatorEnsemble = MRSTHistoryMatchingEnsemble(baseExample, operatorSamples, qoi, ...
    'directory', fullfile(topDirectory, 'operator'), ...
    'simulationStrategy', 'parallel', ...
    'maxWorkers', 8, ...
    'reset', true, ...
    'verbose', true)

%% Run ensemble
operatorEnsemble.simulateEnsembleMembers();


%% Get simulated observations
disp('simulated observations: ')
size(operatorEnsemble.getEnsembleQoI())

%% Get the matrix of ensemble samples 
disp('Matrix of ensemble samples (parameters):')
size(operatorEnsemble.getEnsembleSamples())


%% Do history matching
operatorEnsemble.doHistoryMatching()

%% Run new ensemble
operatorEnsemble.simulateEnsembleMembers();

%% Plot original and updated ensemble results
operatorEnsemble.plotQoI('subplots', true, 'clearFigure', false, ...
    'cmapName', 'lines', ...
    'plotTruth', true, ...
    'legend', {'observations', 'truth', 'posterior mean', 'prior mean'});
