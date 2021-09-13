%% History matching of reduced network model of the Egg model

mrstModule add ad-core ad-blackoil mrst-gui ad-props ...
    example-suite incomp ensemble network-models diagnostics

mrstVerbose off

%% Set the name of the base problem and ensemble size
% The base problem contains all properties that are common throughout the
% ensemble


egg_realization = 99;
baseProblemName = 'egg_wo_network';
trueProblemName = 'egg_wo';

baseProblemOptions = {'realization', egg_realization, ...
                      'fullSchedule', false};

ensembleSize = 100;
%ensembleSize = 160;

%% Define where to store truth and ensemble simulations
directoryTruth = fullfile(mrstOutputDirectory(), ...
                          'historyMatching', 'truth','egg_99', ...
                          baseProblemName);
                      
topDirectory = fullfile(mrstOutputDirectory(), ...
                        'historyMatching', 'tutorial', 'egg_network_99', baseProblemName);

%% Create example
baseExample = MRSTExample(baseProblemName, ...
                          'deleteOldResults', false, ...
                          'plotNetwork', false, ...
                          baseProblemOptions{:});

numConnections = numel(baseExample.options.connectionIndices.cells);

                      
%% Create the full model so that we can generate observations
trueExample = MRSTExample(baseExample.options.fullExampleName);
trueProblem = trueExample.getPackedSimulationProblem('Directory', directoryTruth);

plotExample = false;
rerunTrueProblemFromScratch = false;
overwriteObservation = true;

if rerunTrueProblemFromScratch
    clearPackedSimulatorOutput(trueProblem, 'prompt', false);
end
simulatePackedProblem(trueProblem);
if plotExample
    [wellSols, states, reports] = getPackedSimulatorOutput(trueProblem);
    trueExample.plot(states);
end

relevantWellNames = {baseExample.schedule.control.W.name};
%productionWellNames = productionWellNames(9:end);
numWells = numel(relevantWellNames);

% Read and store the observations in a QoI object.
% We are interested in the fluxes of oil and water from the producers, and
% the bottom hole pressure from the injectors.
trueQoI = WellQoIHM('wellNames', relevantWellNames, ...
                    'names', {'qOs', 'qWs', 'bhp'});
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
obsStdDevFlux = 50*stb()/day();
obsStdDevBhp  = 1*barsa();
obsStdDev = [obsStdDevFlux, obsStdDevFlux, obsStdDevBhp];
trueObservations  = trueQoI.getQoI(trueProblem); 
if numel(observationResultHandler.getValidIds) < 1 || overwriteObservation
    for w = 1:numel(trueQoI.wellNames)
        perturbedObservations(w) = trueObservations(w);
        for f = 1:numel(trueQoI.names)
            trueVals = trueObservations(w).(trueQoI.names{f});
            perturbedObservations(w).(trueQoI.names{f}) = trueVals; %  + randn(size(trueVals))*obsStdDev;
        end
    end
    observationResultHandler{1} = {perturbedObservations};
end
if numel(truthResultHandler.getValidIds) < 1 || overwriteObservation
    truthResultHandler{1} = {trueObservations};
end


                      
%% Define samples


% Uncertain initial conditions (saturation of oil)
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

% Uncertain transmissibility and porevolume
operatorData = cell(ensembleSize, 1);
FDmean = baseExample.model.operators.pv(baseExample.options.connectionIndices.cells{1}(1));
logMean = log(baseExample.model.operators.T(baseExample.options.connectionIndices.faces{1}(1)));
for i = 1:ensembleSize
    operatorData{i}.pv = max((baseExample.options.pv + (baseExample.options.pv/6).*randn(size(baseExample.options.pv))), baseExample.options.pv/3)';
    operatorData{i}.T = exp(log(baseExample.options.TT)+ 0.3*randn(size(baseExample.options.TT)))';
    
    
    %operatorData{i}.pv = max(FDmean + (FDmean/6)*randn(1, numConnections), FDmean/3);
    %operatorData{i}.T = exp(logMean + 1*randn(1, numConnections));
    
    operatorData{i}.pv = operatorData{i}.pv./10;
    %operatorData{i}.T  = operatorData{i}.T.*10;
end


% Uncertain well production indices
wellSampleData = cell(ensembleSize, 1);
for i = 1:ensembleSize
    % wellSampleData{i}.WI = (rand(1,numWells)+0.1)*1e-11;
    % wellSampleData{i}.WI = wellSampleData{i}.WI.*2;
    wellSampleData{i}.WI = (randn(numWells,1)*0.25 + 1.25)*1e-11;   
end

wellSamples = WellSamplesHM('data', wellSampleData); %, ...
                            %'WIScale', 1e-11);

operatorSamples = NetworkOperatorSamplesHM('data', operatorData, ...
                                     'connectionIndices', baseExample.options.connectionIndices, ...
                                     'pvScale', 1e4); %, ...
                                     %'TScale', 1e-9);
                                     %  'pvScale', baseExample.model.G.cells.volumes(1)/100)
                                     

initSoSamples = NetworkState0SamplesHM('data', initSoData, ...
                                       'connectionIndices', baseExample.options.connectionIndices)

%% Making a composite sample object:

%samples = CompositeSamplesHM({wellSamples, operatorSamples, initSoSamples})

samples = CompositeSamplesHM({operatorSamples, wellSamples, initSoSamples});

%samples = CompositeSamplesHM({operatorSamples, wellSamples});



%% Create QoI, with the observations
qoi = WellQoIHM(...
    'wellNames', relevantWellNames, ...
    'names', {'qOs', 'qWs', 'bhp'}, ...
    'observationResultHandler', observationResultHandler, ...
    'truthResultHandler', truthResultHandler, ...
    'observationCov', obsStdDev.^2);




%% Create the ensemble
ensemble = MRSTHistoryMatchingEnsemble(baseExample, samples, qoi, ...
    'alpha', [28/3 7 4 2], ...
    'directory', fullfile(topDirectory, 'composite'), ...
    'simulationStrategy', 'spmd', ...
    'maxWorkers', 8, ...
    'reset', true, ...
    'verbose', true);

totalNumberOfTimesteps = numel(ensemble.originalSchedule.step.val);
ensemble.updateHistoryMatchingInterval(2:2:totalNumberOfTimesteps);

%% Run ensemble
ensemble.simulateEnsembleMembers();

%ensemble.plotQoI('plotTruth', true);


%% Get simulated observations
disp('simulated observations: ')
size(ensemble.getEnsembleQoI())

%% Get the matrix of ensemble samples 
disp('Matrix of ensemble samples (parameters):')
size(ensemble.getEnsembleSamples())


%% Do history matching
ensemble.doHistoryMatching()

%% Run new ensemble at full length
ensemble.updateHistoryMatchingInterval(1:totalNumberOfTimesteps);
ensemble.simulateEnsembleMembers();

%% Plot original and updated ensemble results


ensemble.plotQoI('subplots', false, 'clearFigure', false, ...
    'cmapName', 'lines', ...
    'plotTruth', true, ...
    'subIterations', true, ...
    'observationIndices', (2:2:totalNumberOfTimesteps), ...
    'legend', {'observations', 'truth', 'posterior mean', 'ES-MDA it 3',...
               'ES-MDA it 2', 'ES-MDA it 1', 'prior mean'});

%% Copyright Notice
%
% <html>
% <p><font size="-1">
% Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.
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
