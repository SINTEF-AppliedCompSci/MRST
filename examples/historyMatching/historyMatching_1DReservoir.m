%% History mathing of 1D reservoir simulation
% In this example, we use well production data to estimate 
mrstModule add ad-core ad-blackoil mrst-gui ad-props ...
    example-suite incomp ensemble

mrstVerbose off


%% Set up and simulate the true solution
% We will here use an identical twin experiment, where we use the same
% problem for both generating the truth and as a base for our ensemble.


trueProblemName = 'ensemble_base_problem_1d_reservoir';
numCells = 10;
trueProblemOptions = {'ncells', numCells, ...
                      'rngseed', 1};
                  
directoryTruth = fullfile(mrstOutputDirectory(), ...
                          'historyMatching', 'truth', ...
                          trueProblemName);
                  
trueExample = MRSTExample(trueProblemName, trueProblemOptions{:});
trueProblem = trueExample.getPackedSimulationProblem('Directory', directoryTruth);

plotExample = false;
rerunTrueProblemFromScratch = false;


if rerunTrueProblemFromScratch
    clearPackedSimulatorOutput(trueProblem);
end
simulatePackedProblem(trueProblem);
if plotExample
    [wellSols, states, reports] = getPackedSimulatorOutput(trueProblem);
    trueExample.plot(states);
end



%% Generate observations 
% Define a QOI for storing the relevant qoi for our problem

trueQoI = WellQoI('wellNames', {'P1'}, 'fldname', {'qOs'});
trueQoI = trueQoI.validateQoI(trueProblem);
trueObservations = trueQoI.getQoI(trueProblem);

% Define observation uncertainty 
obsStdDev = 0.0004;

% Create a separate ResultHandler for the observations 
observationResultHandler = trueQoI.ResultHandler;
observationResultHandler.dataPrefix = 'observedQoI';

% Add some observation noise and store output
if numel(observationResultHandler.getValidIds) < 1
    observations = trueObservations{1} + randn(size(trueObservations{1}))*obsStdDev;
    observationResultHandler{1} = {observations};
end

%% Select and populate samples for the stochastic components in the ensemble

ensembleSize = 23;

configData = cell(ensembleSize, 1);
for i = 1:ensembleSize
    configData{i}.poro = gaussianField(trueExample.model.G.cartDims, [0.2 0.4]); 
    configData{i}.perm = configData{i}.poro.^3.*(1e-5)^2./(0.81*72*(1-configData{i}.poro).^2);
end

samples = RockSamples('data', configData);

%% Select quantity of interest class matching the what we have as observations
% We validate the QoI with the trueProblem, since this will be our ensemble
% base problem as well.

qoi = WellQoI('wellNames', {'P1'}, 'fldname', {'qOs'}, ...
                  'observationResultHandler', observationResultHandler, ...
                  'observationCov', obsStdDev^2);

%% Create the ensemble

ensemble = MRSTHistoryMatchingEnsemble(trueExample, samples, qoi, ... 
    ... %'directory', uniqueDirectory, ...
    'simulationType', 'parallel', ...
    'maxWorkers', 8, ...
    'verbose', true, ...
    'deleteOldResults', true...
    );


%% Displaying the observations and observation error cov through the ensemble
disp('observation and scaling vector')
[obsVector, obsScaling] = ensemble.qoi.getObservationAndScaling()
disp('observation error covariance matrix')
ensemble.qoi.getObservationErrorCov()


%% Run ensemble
ensemble.simulateAllEnsembleMembers();

%% Get simulated observations
disp('simulated observations')
size(ensemble.getEnsembleQoI())

%% Get the matrix of ensemble samples 
size(ensemble.getEnsembleSamples())

%% Do history matching
disp('updated sample object:')
updatedSamples = ensemble.doHistoryMatching()


%% Create a new ensemble with updated samples
updatedEnsemble = MRSTHistoryMatchingEnsemble(trueExample, updatedSamples, qoi, ... 
    ... %'directory', uniqueDirectory, ...
    'simulationType', 'parallel', ...
    'maxWorkers', 8, ...
    'verbose', true, ...
    'deleteOldResults', false, ...
    'historyMatchingIteration', 2 ...
    );

%% Run new ensemble
updatedEnsemble.simulateAllEnsembleMembers();

%% Plot first ensemble
ensemble.qoi.plotEnsemble(ensemble);

%% Plot updated ensemble
updatedEnsemble.qoi.plotEnsemble(updatedEnsemble);


%% Plot diff between the qois
figure
hold on
for i = 1:ensemble.num
    plot(updatedEnsemble.qoi.ResultHandler{i}{1} - ensemble.qoi.ResultHandler{i}{1})
end







