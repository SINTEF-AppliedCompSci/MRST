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
    example-suite incomp ensemble network-models diagnostics linearsolvers

mrstVerbose off

%% Set the name of the base problem and ensemble size
% The base problem contains all properties that are common throughout the
% ensemble

baseProblemName = 'network_model_template';
trueProblemName = 'olympus_field_wo';

ensembleSize = 100;
numConnections = 77; % Not a parameter, but hard-coded in the baseProblem

%% Define where to store truth and ensemble simulations
directoryTruth = fullfile(mrstOutputDirectory(), ...
                          'historyMatching', 'truth','network_operator', ...
                          baseProblemName);
                      
topDirectory = fullfile(mrstOutputDirectory(), ...
                        'historyMatching', 'tutorial', 'network_operator', baseProblemName);

%% Create example
baseExample = MRSTExample(baseProblemName, ...
                          'fullExampleName',trueProblemName,...
                          'gpsnetPerm', 500*milli*darcy,...
                          'deleteOldResults', false, ...
                          'plotNetwork', false,...
                          'realization',1);

                      
%% Create the full model so that we can generate observations
trueExample = MRSTExample(trueProblemName);
trueProblem = trueExample.getPackedSimulationProblem('Directory', directoryTruth);
%simulatePackedProblem(trueProblem);

plotExample = true;
rerunTrueProblemFromScratch = false;
overwriteObservation = true;

if rerunTrueProblemFromScratch
    clearPackedSimulatorOutput(trueProblem, 'prompt', true);
end
simulatePackedProblem(trueProblem);
if plotExample
    [wellSols, states, reports] = getPackedSimulatorOutput(trueProblem);
    trueExample.plot(states);
end

 wellNames={baseExample.schedule.control.W.name};%{'INJ-1','INJ-2','INJ-3','INJ-4','INJ-5','INJ-6','INJ-7','PROD-1','PROD-10','PROD-11','PROD-2','PROD-3','PROD-4','PROD-5','PROD-6','PROD-7','PROD-8','PROD-9'};
 numWells = numel(wellNames);

% Read and store the observations in a QoI object
trueQoI = WellQoIHM('wellNames', wellNames, ...
                    'names', {'qOs', 'qWs'});
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
obsStdDev = 0.0004;
trueObservations  = trueQoI.getQoI(trueProblem); 
if numel(observationResultHandler.getValidIds) < 1 || overwriteObservation
    for w = 1:numel(trueQoI.wellNames)
        perturbedObservations(w) = trueObservations(w);
        for f = 1:numel(trueQoI.names)
            trueVals = trueObservations(w).(trueQoI.names{f});
            perturbedObservations(w).(trueQoI.names{f}) = trueVals + 0*randn(size(trueVals))*obsStdDev; % no perturbation
        end
    end
    observationResultHandler{1} = {perturbedObservations};
end
if numel(truthResultHandler.getValidIds) < 1 || overwriteObservation
    truthResultHandler{1} = {trueObservations};
end


                      
%% Define samples that give different transmissibilise for each connection
rng(12345)
% Initializing as a log-Gaussian distribution around the value obtained
% from flow diagnostics 
% (assuming FD finds the same T for all connections)
logMean = log(baseExample.model.operators.T(baseExample.options.connectionIndices.faces{1}(1)));
transData = cell(ensembleSize, 1);
for i = 1:ensembleSize
    transData{i}.T = exp(logMean + 1*randn(1, numConnections));
end


%% Create sample object
transSamples = NetworkOperatorSamplesHM('data', transData,...
                                        'TScale',max(transData{i}.T),...
                                        'connectionIndices', baseExample.options.connectionIndices)
%% Create QoI, with the observations







%% History match with porevolumes
FDmean = 0.3*baseExample.model.operators.pv(baseExample.options.connectionIndices.cells{1}(1));
pvData = cell(ensembleSize, 1);
for i = 1:ensembleSize
    pvData{i}.pv = FDmean + 0.1*FDmean*randn(1, numConnections);
end

%baseExample.schedule.control.W
wellSampleData = cell(ensembleSize, 1);
for i = 1:ensembleSize
    wellSampleData{i}.WI = (randn(numWells,1)*0.5 + 10)*1e-11;   
    %wellSampleData{i}.WI = precomputedSamples.WI_sum(:, i); 
end
minWI = 0.01*min(min(wellSampleData{i}.WI));
maxWI = 50*max(max(wellSampleData{i}.WI));

wellSamples = WellSamplesHM('data', wellSampleData, ...
                            'WIScale', 1e-11, ...
                            'minWIValue', minWI, ...
                            'maxWIValue', maxWI);



%% History match with both porevolumes and transmissibilities
operatorData = cell(ensembleSize, 1);
for i = 1:ensembleSize
    operatorData{i}.pv = pvData{i}.pv;
    operatorData{i}.T = transData{i}.T;
end

operatorSamples = NetworkOperatorSamplesHM('data', operatorData, ...
                                     'connectionIndices', baseExample.options.connectionIndices, ...
                                     'pvScale', baseExample.model.G.cells.volumes(1),...
                                     'TScale',max(transData{i}.T))

wellSamples = WellSamplesHM('data', wellSampleData, ...
                            'WIScale', 1e-11, ...
                            'minWIValue', minWI, ...
                            'maxWIValue', maxWI);
                        
                        
%% Define samples that give different initial saturation for each connection

initSoData = cell(ensembleSize, 1);
for i = 1:ensembleSize
%     initSoTmp = zeros(1,numConnections);
%    while((min(initSoTmp(:)) < 0.8) || (max(initSoTmp(:)) > 1))
        initSoTmp = 0.3*rand(1,numConnections)+0.7;%Uniform distribution
%          initSoTmp = randn(1,numConnections);
%         initSoTmp = initSoTmp - mean(initSoTmp(:));
%         initSoMin = min(initSoTmp(:));
%         initSoMax = max(initSoTmp(:));
%         initSoTmp = initSoTmp.*(0.45/(initSoMax - initSoMin)) + 0.9;
%    end
    
    initSoData{i}.initSo = initSoTmp;
end


%% Create sample object
initSoSample = NetworkState0SamplesHM('data', initSoData, ...
                                      'connectionIndices', baseExample.options.connectionIndices)
                        
%% Making a composite sample object:
samples = CompositeSamplesHM({initSoSample,operatorSamples,wellSamples});
   
%% 
qoi = WellQoIHM(...
    'wellNames', wellNames, ...
    'names', {'qOs', 'qWs'}, ...
    'observationResultHandler', observationResultHandler, ...
    'truthResultHandler', truthResultHandler, ...
    'observationCov', obsStdDev^2);


%% Create combo ensemble
operatorEnsemble = MRSTHistoryMatchingEnsemble(baseExample, samples, qoi, ...
    'alpha', [28/3 7 4 2], ...
    'directory', fullfile(topDirectory, 'operator'), ...
    'simulationStrategy', 'spmd', ...
    'maxWorkers', 20, ...
    'verboseSimulation',false,...
    'reset', true, ...
    'storeOutput', true,...
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
operatorEnsemble.plotQoI('subplots', false, 'clearFigure', false, ...
    'cmapName', 'lines', ...
    'plotTruth', true, ...
    'subIterations', true, ...
    'legend', {'observations', 'truth', 'posterior mean', 'ES-MDA it 3',...
               'ES-MDA it 2', 'ES-MDA it 1', 'prior mean'});


% operatorEnsemble.plotQoI('subplots', false, 'clearFigure', false, ...
%     'cmapName', 'lines', ...
%     'plotTruth', true, ...
%     'legend', {'observations', 'truth', 'prior mean'});
