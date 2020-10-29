%% History matching of 3D reservoir simulation
% This is a minimal example (and a sanity test) for creating an ensemble of
% 3D reservoirs with two phases (oil and water).
mrstModule add ad-core ad-blackoil mrst-gui ad-props ...
    example-suite incomp ensemble 

mrstVerbose off

%% Set up and simulate the true solution
% We will here use an identical twin experiment, where we use the same
% problem for both generating the truth and as a base for our ensemble.

trueProblemName = 'ensemble_base_problem_3d_reservoir';
baseProblemOptions = {};


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
% Define a QoI object for storing the relevant observations we will use for
% history matching

trueQoI = WellQoI(...
    'wellNames', {'P1', 'P2'}, ...
    'fldname', {'qOs', 'qWs'}, ...
    'cumulative', false, ...
    'numTimesteps', []);

trueQoI = trueQoI.validateQoI(trueProblem);
trueObservations = trueQoI.getQoI(trueProblem);

% Define observation uncertainty 
obsStdDev = 0.0004*0.1;

% Create a separate ResultHandler for the observations 
observationResultHandler = trueQoI.ResultHandler;
observationResultHandler.dataPrefix = 'observedQoI';

% Add some observation noise and store output
if numel(observationResultHandler.getValidIds) < 1
    observations = trueObservations{1} + randn(size(trueObservations{1}))*obsStdDev;
    observationResultHandler{1} = {observations};
end



%% Select and populate samples for stochastic configurations class

ensembleSize = 23;


rockData = cell(ensembleSize, 1);
for i = 1:ensembleSize
    rockData{i}.poro = gaussianField(trueExample.model.G.cartDims, [0.2 0.4]); 
    rockData{i}.perm = rockData{i}.poro.^3.*(1e-5)^2./(0.81*72*(1-rockData{i}.poro).^2);
end

rockSamples = RockSamples('data', rockData);

%% Select quantity of interest class matching the what we have as observations
% We validate the QoI with the trueProblem, since this will be our ensemble
% base problem as well.

qoi = WellQoI('wellNames', {'P1', 'P2'}, ...
              'fldname', {'qOs', 'qWs'}, ...
              'observationResultHandler', observationResultHandler, ...
              'observationCov', obsStdDev^2);


%% Create the ensemble
originalRockEnsemble = MRSTHistoryMatchingEnsemble(trueExample, rockSamples, qoi, ...
    ... %'directory', uniqueDirectory, ...
    'simulationType', 'parallel', ...
    'maxWorkers', 8, ...
    'deleteOldResults', true, ...
    'verbose', true);

%% Displaying the observations and observation error cov through the ensemble
disp('observation and scaling vector')
[obsVector, obsScaling] = originalRockEnsemble.qoi.getObservationAndScaling()
disp('observation error covariance matrix')
originalRockEnsemble.qoi.getObservationErrorCov()

%% Run ensemble
originalRockEnsemble.simulateAllEnsembleMembers();

%% Get simulated observations
disp('simulated observations')
size(originalRockEnsemble.getEnsembleQoI())

%% Get the matrix of ensemble samples 
size(originalRockEnsemble.getEnsembleSamples())

%% Do history matching
disp('updated sample object:')
updatedRockSamples = originalRockEnsemble.doHistoryMatching()



%% Create a new ensemble with updated samples
updatedRockEnsemble = MRSTHistoryMatchingEnsemble(trueExample, updatedRockSamples, qoi, ...
    ... %'directory', uniqueDirectory, ...
    'simulationType', 'parallel', ...
    'maxWorkers', 8, ...
    'historyMatchingIteration', 2 ...
    );



%% Run new ensemble
updatedRockEnsemble.simulateAllEnsembleMembers();

%% Plot original and updated ensemble results
originalRockEnsemble.qoi.plotEnsemble(originalRockEnsemble);

updatedRockEnsemble.qoi.plotEnsemble(updatedRockEnsemble);

%% Investigate change in samples 

originalRockParameters = originalRockEnsemble.getEnsembleSamples();
updatedRockParameters = updatedRockEnsemble.getEnsembleSamples();

originalWellQoI = originalRockEnsemble.getEnsembleQoI();
updatedWellQoI = updatedRockEnsemble.getEnsembleQoI();


%% 

figure
hold on
for i = 1:ensemble.num
    plot(updatedRockParameters(:,i) - originalRockParameters(:,i))
end

figure
hold on
for i = 1:ensemble.num
    plot(updatedWellQoI(:,i) - originalWellQoI(:,i))
end










%% Create another ensemble using stochastic well indices
% ---------------------
wellSampleData = cell(ensembleSize, 1);
for i = 1:ensembleSize
    wellSampleData{i}.WI = rand(1,4)*1e-11;
end

wellSamples = WellSamples('data', wellSampleData);

%% Define new ensemble
originalWellEnsemble = MRSTHistoryMatchingEnsemble(trueExample, wellSamples, qoi, ...
    ... %'directory', uniqueDirectory, ...
    'simulationType', 'parallel', ...
    'maxWorkers', 8, ...
    'deleteOldResults', true, ...
    'verbose', true);

%% Simulate
originalWellEnsemble.simulateAllEnsembleMembers();

%% Do history matching
updatedWellSamples = originalWellEnsemble.doHistoryMatching()

%% Create a new ensemble with updated samples
updatedWellEnsemble = MRSTHistoryMatchingEnsemble(trueExample, updatedWellSamples, qoi, ...
    ... %'directory', uniqueDirectory, ...
    'simulationType', 'parallel', ...
    'maxWorkers', 8, ...
    'historyMatchingIteration', 2 ...
    );

%% Run new ensemble
updatedWellEnsemble.simulateAllEnsembleMembers();

%% Plot original and updated ensemble results
originalWellEnsemble.qoi.plotEnsemble(originalWellEnsemble);

updatedWellEnsemble.qoi.plotEnsemble(updatedWellEnsemble);

%% Investigate change in samples 

originalWellParameters = originalWellEnsemble.getEnsembleSamples();
updatedWellParameters = updatedWellEnsemble.getEnsembleSamples();

originalWellQoI = originalWellEnsemble.getEnsembleQoI();
updatedWellQoI = updatedWellEnsemble.getEnsembleQoI();



figure
hold on
for i = 1:ensemble.num
    plot(updatedWellParameters(:,i) - originalWellParameters(:,i))
end

figure
hold on
for i = 1:ensemble.num
    plot(updatedWellQoI(:,i) - originalWellQoI(:,i))
end












%% Create an ensemble that combines both sampling strategies
% ---------------------

% We first organize our precomputed samples in a struct that is supported
% by the combined sample class
comboData = cell(ensembleSize, 1);
for i = 1:ensembleSize
    comboData{i}.rock = rockData{i};
    comboData{i}.well = wellSampleData{i};
end

comboSamples = WellRockSamples('data', comboData);

%% Define new ensemble
originalComboEnsemble = MRSTHistoryMatchingEnsemble(trueExample, comboSamples, qoi, ...
    ... %'directory', uniqueDirectory, ...
    'simulationType', 'parallel', ...
    'maxWorkers', 8, ...
    'deleteOldResults', true, ...
    'verbose', true);

%% Simulate and plot
originalComboEnsemble.simulateAllEnsembleMembers();

%% Do history matching
updatedComboSamples = originalComboEnsemble.doHistoryMatching()

%% Create a new ensemble with updated samples
updatedComboEnsemble = MRSTHistoryMatchingEnsemble(trueExample, updatedComboSamples, qoi, ...
    ... %'directory', uniqueDirectory, ...
    'simulationType', 'parallel', ...
    'maxWorkers', 8, ...
    'historyMatchingIteration', 2 ...
    );

%% Run new ensemble
updatedComboEnsemble.simulateAllEnsembleMembers();

%% Plot original and updated ensemble results
originalComboEnsemble.qoi.plotEnsemble(originalComboEnsemble);

updatedComboEnsemble.qoi.plotEnsemble(updatedComboEnsemble);
