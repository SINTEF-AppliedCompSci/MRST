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
                      
topDirectory = fullfile(mrstOutputDirectory(), ...
                        'historyMatching', 'esmdaTutorial', trueProblemName);
                      
trueExample = MRSTExample(trueProblemName);
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
    'cumulative', false);

trueQoI = trueQoI.validateQoI(trueProblem);
trueObservations = trueQoI.getQoI(trueProblem);

% Define observation uncertainty 
obsStdDev = 0.0004; %*0.1;

% Create a separate ResultHandler for the observations 
observationResultHandler = trueQoI.ResultHandler;
observationResultHandler.dataPrefix = 'observedQoI';

% Add some observation noise and store output
if numel(observationResultHandler.getValidIds) < 1
    for w = 1:numel(trueQoI.wellNames)
        for f = 1:numel(trueQoI.fldname)
            perturbedObservations{w}{f} = trueObservations{w}{f} + randn(size(trueObservations{w}{f}))*obsStdDev;
        end
    end
    observationResultHandler{1} = {perturbedObservations};
end



%% Select and populate samples for stochastic configurations class

ensembleSize = 70;


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
esmdaEnsemble = MRSTHistoryMatchingEnsemble(trueExample, rockSamples, qoi, ...
    'alpha', [28/3 7 4 2], ...
    'directory', fullfile(topDirectory, 'esmda'), ...
    'simulationStrategy', 'parallel', ...
    'maxWorkers', 8, ...
    'reset', true, ...
    'verbose', true)

%% Displaying the observations and observation error cov through the ensemble
disp('observation and scaling vector')
[obsVector, obsScaling] = esmdaEnsemble.qoi.getObservationAndScaling()
disp('observation error covariance matrix')
esmdaEnsemble.qoi.getObservationErrorCov()

%% Run ensemble
esmdaEnsemble.simulateEnsembleMembers();

%% Get simulated observations
disp('simulated observations')
size(esmdaEnsemble.getEnsembleQoI())

%% Get the matrix of ensemble samples 
size(esmdaEnsemble.getEnsembleSamples())

%% Do history matching
disp('updated sample object:')
esmdaEnsemble.doHistoryMatching();

%% Run resulting forecast
esmdaEnsemble.simulateEnsembleMembers();


%% Plot results
esmdaEnsemble.plotQoI('subplots', true, ...
                      'clearFigure', false, ...
                      'subIterations', true, ...
                      'legend', {'observations', 'posterior', 'ES-MDA it 3',...
                                 'ES-MDA it 2', 'ES-MDA it 1', 'prior'});







%% EnKF
%%%%%%%%%%%%%%%%%%%%%

enkfEnsemble = MRSTHistoryMatchingEnsemble(trueExample, rockSamples, qoi, ...
    'alpha', 1, ... %'alpha', [2 2] ,... %[28/3 7 4 2], ...
    'directory', fullfile(topDirectory, 'enkf'), ...
    'simulationStrategy', 'parallel', ...
    'maxWorkers', 8, ...
    'reset', true, ...
    'verbose', true, ...
    'verboseSimulation', false)


%enkfEnsemble.simulateEnsembleMembers('range', 1);


%% Loop 

totalNumTimesteps = 15; %numel(trueProblem.SimulatorSetup.schedule.step.val);
enkfInvervals = {(1:2), (3:4), (5:6), (7:8), (9:10)};

for i = 1:numel(enkfInvervals)
    enkfEnsemble.updateHistoryMatchingInterval(enkfInvervals{i});
    enkfEnsemble.simulateEnsembleMembers();
    enkfEnsemble.doHistoryMatching();
end

%% Forecast the posterior
enkfEnsemble.updateHistoryMatchingInterval(1:15);
enkfEnsemble.simulateEnsembleMembers();

%% Plot history matching

enkfEnsemble.plotQoI('subplots', true, ...
                     'clearFigure', false, ...
                     'subIterations', true, ...
                     'observationIndices', (1:10) ); %, ...
                     
                      





