%% History matching of 3D reservoir simulation
% This is a minimal example (and a sanity test) for creating an ensemble of
% 3D reservoirs with two phases (oil and water).
mrstModule add ad-core ad-blackoil mrst-gui ad-props ...
    test-suite incomp ensemble 

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
                        'historyMatching', 'tutorial', trueProblemName);
                      
baseCase = TestCase(trueProblemName);
baseProblem = baseCase.getPackedSimulationProblem('Directory', directoryTruth);

plotExample = false;
rerunTrueProblemFromScratch = true;
overwriteObservation = true;


if rerunTrueProblemFromScratch
    clearPackedSimulatorOutput(baseProblem, 'prompt', false);
end
simulatePackedProblem(baseProblem);
if plotExample
    [wellSols, states, reports] = getPackedSimulatorOutput(baseProblem);
    baseCase.plot(states);
end



%% Select and populate samples for stochastic configurations class

% Small ensemble given as default for faster run-time
%ensembleSize = 70;
ensembleSize = 8;


rockData = cell(ensembleSize, 1);
for i = 1:ensembleSize
    rockData{i}.poro = gaussianField(baseCase.model.G.cartDims, [0.2 0.4]); 
    rockData{i}.perm = rockData{i}.poro.^3.*(1e-5)^2./(0.81*72*(1-rockData{i}.poro).^2);
end

rockSamples = RockSamplesHM('data', rockData);

%% Select quantity of interest class
% We validate the QoI with the trueProblem, since this will be our ensemble
% base problem as well.

obsStdDev = 0.0004*0.1;

qoi = WellQoIHM('wellNames', {'P1', 'P2'}, ...
                'names', {'qOs', 'qWs'}, ...
                'observationCov', obsStdDev^2);
                %'observationResultHandler', observationResultHandler, ...
                %'truthResultHandler', truthResultHandler, ...
                

%% Create the ensemble
% and specify from where we want to take observations
rockEnsemble = MRSTHistoryMatchingEnsemble(baseCase, rockSamples, qoi, ...
    'observationProblem', baseProblem, ...
    'directory', fullfile(topDirectory, 'rock'), ...
    'simulationStrategy', 'background', ...
    'maxWorkers', 4, ...
    'reset', true, ...
    'verbose', true)

%% Displaying the observations and observation error cov through the ensemble
disp('observation and scaling vector')
[obsVector, obsScaling] = rockEnsemble.qoi.getObservationAndScaling()
disp('observation error covariance matrix')
rockEnsemble.qoi.getObservationErrorCov()

%% Run ensemble
rockEnsemble.simulateEnsembleMembers();

%% Get simulated observations
disp('simulated observations')
size(rockEnsemble.getEnsembleQoI())

%% Get the matrix of ensemble samples 
size(rockEnsemble.getEnsembleSamples())

%% Do history matching and thereby update the samples in the ensemble
disp('updated sample object:')
rockEnsemble.doHistoryMatching()


%% Run ensemble with new samples
rockEnsemble.simulateEnsembleMembers();

%% Plot original and updated ensemble results
%h = rockEnsemble.plotQoI('color', [0 0 1], 'subplots', true);

rockEnsemble.plotQoI('subplots', true, 'clearFigure', false, ...
    'plotTruth', true, ...
    'legend', {'observations', 'truth', 'posterior mean', 'prior mean'});








%% Create another ensemble using stochastic well indices
% ---------------------
wellSampleData = cell(ensembleSize, 1);
for i = 1:ensembleSize
    wellSampleData{i}.WI = rand(1,4)*1e-11;
end

wellSamples = WellSamplesHM('data', wellSampleData);

%% Define new ensemble
wellEnsemble = MRSTHistoryMatchingEnsemble(baseCase, wellSamples, qoi, ...
    'observationProblem', baseProblem, ...
    'directory', fullfile(topDirectory, 'well'), ...
    'simulationStrategy', 'background', ...
    'maxWorkers', 4, ...
    'reset', true, ...
    'verbose', true);

%% Simulate
wellEnsemble.simulateEnsembleMembers();

%% Do history matching
wellEnsemble.doHistoryMatching()

%% Run new ensemble
wellEnsemble.simulateEnsembleMembers();

%% Plot original and updated ensemble results
wellEnsemble.plotQoI('subplots', true, 'clearFigure', false, ...
    'plotTruth', true, ...
    'legend', {'observations', 'truth', 'posterior mean', 'prior mean'});










%% Create an ensemble that combines both sampling strategies
% ---------------------

comboSamples = CompositeSamplesHM({rockSamples, wellSamples}, 'tensorProduct', false)


%% Define new ensemble
comboEnsemble = MRSTHistoryMatchingEnsemble(baseCase, comboSamples, qoi, ...
    'observationProblem', baseProblem, ...
    'directory', fullfile(topDirectory, 'combo'), ...
    'simulationStrategy', 'background', ...
    'maxWorkers', 4, ...
    'reset', true, ...
    'verbose', true);

%% Simulate and plot
comboEnsemble.simulateEnsembleMembers();

%% Do history matching
comboEnsemble.doHistoryMatching()

%% Run new ensemble
comboEnsemble.simulateEnsembleMembers();

%% Plot original and updated ensemble results
comboEnsemble.plotQoI('subplots', true, 'clearFigure', false, ...
    'plotTruth', true, ...
    'legend', {'observations', 'truth', 'posterior mean', 'prior mean'});

%% Copyright Notice
%
% <html>
% <p><font size="-1">
% Copyright 2009-2022 SINTEF Digital, Mathematics & Cybernetics.
% </font></p>
% <p><font size="-1">
% This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).
% </font></p>
% <p><font size="-1">
% MRST is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% </font></p>
% <p><font size="-1">
% MRST is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% </font></p>
% <p><font size="-1">
% You should have received a copy of the GNU General Public License
% along with MRST.  If not, see
% <a href="http://www.gnu.org/licenses/">http://www.gnu.org/licenses</a>.
% </font></p>
% </html>
