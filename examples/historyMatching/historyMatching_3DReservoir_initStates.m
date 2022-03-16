%% History matching of 3D reservoir simulation
% This is a minimal example (and a sanity test) for creating an ensemble of
% 3D reservoirs with two phases (oil and water), when the uncertain
% parameters are initial saturation of oil defined through the
% State0Samples class
mrstModule add ad-core ad-blackoil mrst-gui ad-props ...
    test-suite incomp ensemble 

mrstVerbose off

%% Set up and simulate the true solution
% We will here use an identical twin experiment, where we use the same
% problem for both generating the truth and as a base for our ensemble.

trueProblemName = 'ensemble_base_problem_3d_reservoir';
baseProblemOptions = {};


directoryTruth = fullfile(mrstOutputDirectory(), ...
                          'historyMatching', 'truth', 'state0', ...
                          trueProblemName);
                      
topDirectory = fullfile(mrstOutputDirectory(), ...
                        'historyMatching', 'tutorial', 'state0', trueProblemName);
                      
trueExample = TestCase(trueProblemName);
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

%% Generate observations
% Define a QoI object for storing the relevant observations we will use for
% history matching

trueQoI = WellQoIHM(...
    'wellNames', {'P1', 'P2'}, ...
    'names', {'qOs', 'qWs'}, ...
    'cumulative', false);

trueQoI = trueQoI.validateQoI(trueProblem);
trueObservations = trueQoI.getQoI(trueProblem);

% Define observation uncertainty 
obsStdDev = 0.0004*0.1;

% Create a separate ResultHandler for the observations 
observationResultHandler = trueQoI.ResultHandler;
observationResultHandler.dataPrefix = 'observedQoI';

% Add some observation noise and store output
if numel(observationResultHandler.getValidIds) < 1 || overwriteObservation
    for w = 1:numel(trueQoI.wellNames)
        perturbedObservations(w) = trueObservations(w);
        for f = 1:numel(trueQoI.names)
            trueVals = trueObservations(w).(trueQoI.names{f});
            perturbedObservations(w).(trueQoI.names{f}) = trueVals + randn(size(trueVals))*obsStdDev;
        end
    end
    observationResultHandler{1} = {perturbedObservations};
end



%% Select and populate samples for stochastic configurations class

% Small ensemble given as default for faster run-time
%ensembleSize = 70;
ensembleSize = 8;

% We model the transmisibility as a log-Gaussian distribution with the same
% mean and standard deviation as in trueExample

initSoData = cell(ensembleSize, 1);
for i = 1:ensembleSize
    %initSoTmp = zeros(trueExample.model.G.cartDims);
    %while((min(initSoTmp(:)) < 0.2) || (max(initSoTmp(:)) > 1))
        %initSoTmp = randn(trueExample.model.G.cartDims);
        %initSoTmp = initSoTmp - mean(initSoTmp(:));
        %initSoMin = min(initSoTmp(:));
        %initSoMax = max(initSoTmp(:));
        %initSoTmp = initSoTmp.*(0.45/(initSoMax - initSoMin)) + 0.7;
    %end
    
    %initSoData{i}.initSo = initSoTmp;
    initSoData{i}.initSo = gaussianField(trueExample.model.G.cartDims, [0.5 1]); 
end

initSoSample = State0SamplesHM('data', initSoData)

%% Select quantity of interest class matching the what we have as observations
% We validate the QoI with the trueProblem, since this will be our ensemble
% base problem as well.

qoi = WellQoIHM('wellNames', {'P1', 'P2'}, ...
              'names', {'qOs', 'qWs'}, ...
              'observationResultHandler', observationResultHandler, ...
              'observationCov', obsStdDev^2);


%% Create the ensemble
initSoEnsemble = MRSTHistoryMatchingEnsemble(trueExample, initSoSample, qoi, ...
    'directory', fullfile(topDirectory, 'initSo'), ...
    'simulationStrategy', 'background', ...
    'maxWorkers', 4, ...
    'reset', true, ...
    'verbose', true)

%% Displaying the observations and observation error cov through the ensemble
disp('observation and scaling vector')
[obsVector, obsScaling] = initSoEnsemble.qoi.getObservationAndScaling()
disp('observation error covariance matrix')
size(initSoEnsemble.qoi.getObservationErrorCov())

%% Run ensemble
initSoEnsemble.simulateEnsembleMembers();


%% Get simulated observations
disp('simulated observations')
size(initSoEnsemble.getEnsembleQoI())

%% Get the matrix of ensemble samples 
size(initSoEnsemble.getEnsembleSamples())

%% Do history matching and thereby update the samples in the ensemble
disp('updated sample object:')
initSoEnsemble.doHistoryMatching()


%% Run ensemble with new samples
initSoEnsemble.simulateEnsembleMembers();

%% Plot original and updated ensemble results

initSoEnsemble.plotQoI('subplots', true, 'clearFigure', false, ...
    'legend', {'observations', 'posterior mean', 'prior mean'});

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
