%% Egg: Ensemble-based calibration of upscaled models using ES-MDA
% The aim of this example is to calibrate an ensemble of rich network
% models represented by very coarse grid models obtained through upscaling,
% so that they match reference well data from a full Egg model. To do so,
% we generate coarse models using upscaling on each of the realizations
% from the Egg ensemble. This ensemble will then act as our prior, and we 
% use ES-MDA to obtain a posterior parameter distribution. We consider the 
% mean of these parameter as the calibrated model parameters.
% 
% The model is a water-oil reservoir, and we calibrate pore volumes,
% transmissibilities, and well production indices.

mrstModule add ad-core ad-blackoil mrst-gui ad-props ...
    example-suite incomp ensemble optimization upscaling coarsegrid

mrstVerbose off

warning('This example requires some time to run.')

% Some options for running the example
rerunReferenceModel = true;    
plotReferenceModel = true;
regenerateInitialEnsemble = false;

%% Define reference realization and where to store reference data
% Here, we choose which realization from the Egg ensemble we use as
% reference data, and define where to store the data generated from this
% example.

referenceEggRealization = 0;

topDirectory = fullfile(mrstOutputDirectory(), 'upscaled_network_models', 'eggEnsembleESMDACalibration');
referenceDirectory = fullfile(topDirectory, 'random_schedule_truth', ...
    ['egg_realization_', num2str(referenceEggRealization)]);
fullEnsembleDirectory = fullfile(topDirectory, 'full_models');
historyMatchingDirectory = fullfile(topDirectory, ...
    ['upscaled_ensemble_egg_', num2str(referenceEggRealization)]);
originalScheduleDirectory = fullfile(topDirectory, 'original_schedule_truth', ...
    ['egg_realization_', num2str(referenceEggRealization)]);

%% Setup full 3D reference model
% The reference model is one of the realizations from the Egg ensemble. It
% consists of a 3D oil-water model of a highly chanalized reservoir, driven
% by four producer wells and eight injectors.

referenceExample = MRSTExample('egg_wo', 'realization', referenceEggRealization);

% We consider only the first 48 time steps from the reference model, since
% most of the model dynamics play out during that periode.
numTotalTimesteps = 48;

%referenceExample.schedule = simpleSchedule(referenceExample.schedule.step.val(5:53), ...
%                                           'W', referenceExample.schedule.control.W);
referenceExample.schedule = simpleSchedule(referenceExample.schedule.step.val(1:48), ...
                                           'W', referenceExample.schedule.control.W);
    
%% Random schedule

originalSchedule = referenceExample.schedule;
schedule = referenceExample.schedule;
schedule.step.control= ceil(0.25*(1:numel(schedule.step.val))');

rng(0)
for n=2:max(schedule.step.control)
    schedule.control(n) = schedule.control(1);
    for i=1:numel(schedule.control(n).W)
        W = schedule.control(n).W(i);
        switch W.type
            case 'rate'
                W.val = (.75 + .5*rand)*W.val;
            case 'bhp'
                %if rand < 0.2
                %    W.status = false;
                % else
                    W.val = (.95 + 0.1*rand)*W.val;
                %end
        end
        schedule.control(n).W(i) = W;
    end
end
referenceExample.schedule = schedule;

%% Pack and run reference problem
referenceProblem = referenceExample.getPackedSimulationProblem('Directory', referenceDirectory);

if rerunReferenceModel
    clearPackedSimulatorOutput(referenceProblem, 'prompt',  false);
end
simulatePackedProblem(referenceProblem);

if plotReferenceModel
    [refWellSols, refStates, refReports] = getPackedSimulatorOutput(referenceProblem);
    referenceExample.plot(refStates);
    plotWellSols(refWellSols);
end

%% Organize observations from the well solutions of the reference model
% To do this, we build a QoI (quantity of interest) object as defined in
% the ensemble module of MRST. The relevant data is injection/production
% well rates, as well as BHP of the injectors.
% Since we will use only a subset of the reference data (the
% "observations") for calibration, we store the reference data twice,
% enabling us to plot both the observations and the full reference data
% (the "truth") easily within the ensemble module.

wellNames = {referenceExample.schedule.control(1).W.name};
numWells = numel(wellNames);
fieldNames = {'qOs', 'qWs', 'bhp'};

% Define wells and output variables
observationQoI = WellQoIHM('wellNames', wellNames, ...
                           'names', fieldNames);

% Configure the QoI object to the referenceProblem
observationQoI = observationQoI.validateQoI(referenceProblem);
% Evaluate the QoI from the reference problem
referenceObservations = observationQoI.getQoI(referenceProblem);

% Rename the data prefix for the reference observations, to make it easier
% browse through the data related to this example
observationQoI.ResultHandler.dataPrefix = 'observedQoI';

% Fill the the first resulthandler with the reference data
observationQoI.ResultHandler{1} = {referenceObservations};

% Store the reference data in a dedicated result handler for comparing
% predictive performance of the calibrated models 
truthResultHandler = ResultHandler('dataPrefix', 'trueQoI', ...
                                   'writeToDisk', observationQoI.ResultHandler.writeToDisk,...
                                   'dataDirectory', observationQoI.ResultHandler.dataDirectory, ...
                                   'dataFolder', observationQoI.ResultHandler.dataFolder, ...
                                   'cleardir', false);
truthResultHandler{1} = {referenceObservations};

%% Define the coarse upscaled model
% We create a rich network in the form of a coarse model which we will 
% calibrate according to the reference data. 
%
% The model is structured according to the MRSTExample class, which
% integrates directly into the esemble module.
% For details, see the function upscaled_coarse_network

baseUpscaledModel = MRSTExample('upscaled_coarse_network', ...
                                'partition', [6 6 1], ...
                                'referenceExample', referenceExample, ...
                                'plotCoarseModel', true);
                            
                            
%% Build initial ensemble of coarse models through upscaling
% To generate the prior ensemble, we make a coarse model using upscaling of
% each realization of the Egg ensemble. We then store transmissibilities
% and well production indices. 
%
% Since the porosity of the Egg model is constant for all ensemble members,
% and the coarse partitioning is the same, we get the same pore volumes for
% all coarse realizations. This means that the pore volumes won't be
% updating in the history matching. We therefore add some random
% perturbation to the pore volumes to create some initial spread, and
% thereby allow these parameters to be calibrated as well. 

eggRealizations = [1:100];
ensembleSize = numel(eggRealizations);

pvMin = min(baseUpscaledModel.model.operators.pv)/1000;
transMin = eps;
numInternalFaces = numel(baseUpscaledModel.model.operators.T);
numCells = baseUpscaledModel.model.G.cells.num;
numWells = numel(baseUpscaledModel.schedule.control(1).W);

initParamFolder = fullfile(mrstOutputDirectory(), ...
                         'eggEnsembleCoarseESMDACalibration_examples');
initParamFile = fullfile(initParamFolder, ...
                         ['upscaled_parameters.mat']);

% To save processing time, we store the initial parameters in a .mat file
% so that we don't need to do the upscaling for each time we run this
% example.
if exist(initParamFile, 'file') && ~regenerateInitialEnsemble
    
    fprintf('Read pre-computed initial ensemble from file\n');
    % Read pre-computed initial ensemble
    preComputedEnsemble = load(initParamFile);
    transmissibilities = preComputedEnsemble.transmissibilities;
    poreVolumes = preComputedEnsemble.poreVolumes;
    pertPoreVolumes = preComputedEnsemble.pertPoreVolumes;
    %wellProductionIndices = preComputedEnsemble.wellProductionIndices;
    
    wellProductionIndices = preComputedEnsemble.WImean;
    
    % Read transmissibilities from the old prior    
    %oldParams = load('C:\Users\havardh\Documents\MATLAB\mrst-bitbucket\EggFullUpscaleEnsemble.mat');
    %transmissibilities = oldParams.Transmisibility';
    %poreVolumes = oldParams.PoreVolume';
    %wellProductionIndices = oldParams.WI_mean';
    
else
    % Allocate matrices for temporal storage of ensemble parameters
    poreVolumes = ones(ensembleSize, numCells).*pvMin;
    pertPoreVolumes = ones(ensembleSize, numCells).*pvMin;
    wellProductionIndices = zeros(ensembleSize, numWells);
    transmissibilities = ones(ensembleSize, numInternalFaces).*transMin;
    
    WIsum = zeros(ensembleSize, numWells);
    WImean = zeros(ensembleSize, numWells);
    WIharmonic = zeros(ensembleSize, numWells);
    

    totalVolume = sum(baseUpscaledModel.model.operators.pv);
    fractionsOfTotalVolume = baseUpscaledModel.model.operators.pv/totalVolume;
    
    
    for eggRealization = eggRealizations
        fullRealization = MRSTExample('egg_wo', 'realization', eggRealization);
        
        % Use the same time steps as the reference model, but keep the well
        % configurations
        fullRealization.schedule = simpleSchedule(referenceExample.schedule.step.val, ...
                                                  'W', fullRealization.schedule.control.W);

        %fullRealization.schedule.control.W(1).WI
        
        coarseRealization = MRSTExample('upscaled_coarse_network', ...
                                        'partition', [6 6 1], ...
                                        'referenceExample', fullRealization, ...
                                        'plotCoarseModel', false);

        poreVolumes(eggRealization, :) = coarseRealization.model.operators.pv(:);
        transmissibilities(eggRealization, :) = coarseRealization.model.operators.T(:);
        wellProductionIndices(eggRealization, :) = [coarseRealization.schedule.control.W.WI]';

        % Sample random numbers for perturbing pore volumes      
        pvPerturbations = randn(numCells,1).*sqrt(coarseRealization.model.operators.pv);
        
        % Ensure that the perturbation does not alter the total available
        % reservoir volume
        pertVolume = sum(pvPerturbations);
        pvPerturbations = pvPerturbations - fractionsOfTotalVolume.*pertVolume;
        pertPoreVolumes(eggRealization, :) = pvPerturbations;
        
        % Store the different well index upscaling methods
        WIsum(eggRealization, :) = [coarseRealization.options.sum]';
        WImean(eggRealization, :) = [coarseRealization.options.mean]';
        WIharmonic(eggRealization, :) = [coarseRealization.options.harmonic]';
        
        fprintf('Upscaling realization %d%% complete\n', eggRealization);
    
    end     
    
    % Store parameters:
    if ~exist(initParamFolder, 'dir')
        mkdir(initParamFolder);
    end
    save(initParamFile, 'transmissibilities', 'wellProductionIndices', ...
         'poreVolumes', 'pertPoreVolumes', ...
         'WIsum', 'WImean', 'WIharmonic');
end

%% Perturb pore volumes 
% Add the perturbations to the pore volumes at the proper scale.
pertStd = 10;
poreVolumes = poreVolumes + pertPoreVolumes*pertStd;

%% Create sample object from initial ensemble
% We store the initial ensemble parameters in the appropriate Sample
% objects designed for history matching through the ensemble module. Both
% pore volumes and transmissibilities belong to the operator sample class,
% whereas well productivity index is found in well samples. We also specify 
% maximum and minimum allowed values to avoid that simulations break down 
% due to nonphysical configurations.

% The parameters transmissibility and pore volume belongs to the operators
% of a model.
operatorData = cell(ensembleSize, 1);
for i = 1:ensembleSize
    operatorData{i}.pv = poreVolumes(i, :)';
    operatorData{i}.T = transmissibilities(i, :)';
end
transMax = 10*max(max(transmissibilities));
pvMax = 10*max(max(poreVolumes));
operatorSamples = OperatorSamplesHM('data', operatorData, ...
    ... % parameter scaling 
    'pvScale', 1e4, 'TScale', 1e-9, ...
    'minPvValue', pvMin, 'maxPvValue', pvMax, ...
    'minTValue', transMin, 'maxTValue', transMax);

% Well production indices is considered a well parameter
wellSampleData = cell(ensembleSize, 1);

for i = 1:ensembleSize
    wellSampleData{i}.WI = wellProductionIndices(i, :)*2;
end
minWI = 0.01*min(min(wellProductionIndices(:,:)));
maxWI = 8*max(max(wellProductionIndices(:,:)));

wellSamples = WellSamplesHM('data' ,wellSampleData, ...
                            'WIScale', 1e-11, ...
                            'minWIValue', minWI, 'maxWIValue', maxWI);
                        
% Wrap the two sample objects in a single composite object
samples = CompositeSamplesHM({operatorSamples, wellSamples});

%% Create QoI object and define observation uncertainty
% The quantity of interest for an ensemble simulation defines which data we
% will keep from after simulating each realization, and in a history
% matching context, which data to use for calibration. The QoI here should
% therefore follow the same structure as the reference data that we defined
% in the start of this example.
% 
% In history matching, the observations are always considered to contain
% some noise, or uncertainty, and history matching algorithms, such as
% ES-MDA, is constructed to calibrate model parameters with respect to the
% uncertainty of the observations. In this example, we consider the
% reference data to be exact, but we still need to give a measure of the
% observation uncertainty for the algorithm to make sense. 

% Observation uncertainty
obsStdDevFlux = 50*stb()/day();
obsStdDevBhp  = 1*barsa();
obsStdDev = [obsStdDevFlux, obsStdDevFlux, obsStdDevBhp];

% Quantity of interest
qoi = WellQoIHM('wellNames', wellNames, ...
                'names', fieldNames, ...
                'observationResultHandler', observationQoI.ResultHandler, ...
                'truthResultHandler', truthResultHandler, ...
                'observationCov', obsStdDev.^2);
            
%% Create the ensemble object
% Here, we gather all the information we have made until now. The ensemble
% consists of the baseExample that defines everything all the ensemble
% members have in common (geological model, grid, fluid, simulation
% methods, and so on). We also specify the alpha parameters for ES-MDA,
% choose parallel simulation strategy, etc. More information about input
% parameters can be found in the base class MRSTEnsemble.m.

% The maxWorkers argument determines the maximum number of simultaneous 
% processes which will be used to simulate the ensemble. Increasing this
% number will increase the number of processes used, but care must be taken
% to ensure that your computer is capable of handling this number of 
% processes. Specifying too high a number will result in high memory usage 
% which could drastically slow down your computer.

ensemble = MRSTHistoryMatchingEnsemble(baseUpscaledModel, samples, qoi, ...
    'alpha', [28/3 7 4 2], ...
    'directory', historyMatchingDirectory, ...
    'simulationStrategy', 'spmd', ...
    'maxWorkers', 8, ...
    'reset', true, ...
    'verbose', true, ...
    'verboseSimulation', false);

%% Run prior ensemble 
% If simulations are run in parallel (i.e. maxWorkers > 1) but you do not
% have the MATLAB Parallel toolbox, then maxWorkers number of individual
% MATLAB processes will be started in the background. These processes will
% be killed automatically when the simulation is finished, however, if the
% simulation fails then you will need to manually kill these processes via
% your system task manager.

ensemble.simulateEnsembleMembers('progressTitle', 'Simulating prior ensemble');


%% Configure the schedule used during history matching
% During the calibration, we only use every second observation, and we only
% consider the first half of the simulation time span. The second half will
% be used to look at predictive qualities of the results.

totalNumberOfTimesteps = numel(ensemble.originalSchedule.step.val);
%observationIndices = (2:2:floor(totalNumberOfTimesteps/2));
observationIndices = (1:totalNumberOfTimesteps);
ensemble.updateHistoryMatchingInterval(observationIndices);

%% Do history matching
% This function performs ES-MDA updates with the relaxation parameters
% given as 'alpha' to the ensemble, running intermediate ensemble 
% simulations as it goes.
ensemble.doHistoryMatching()

%% Run the posterior ensemble using the calibrated parameters
% We specify that we want to run the posterior for the complete simulation
% time span.
ensemble.updateHistoryMatchingInterval(1:totalNumberOfTimesteps);
ensemble.simulateEnsembleMembers('progressTitle', 'Simulating posterior ensemble');


%% Plot prior and posterior ensemble results
%
% Create a boolean array of structs that reflects the QoI object, where
% each boolean value reflects which
plotWells = [];
for w = 1:numWells
    for f = 1:numel(ensemble.qoi.names)
        plotWells(w).(ensemble.qoi.names{f}) = false;
    end
end

% Plot bhp from wells 2 and 5, and the production rates from well 12 (PROD4)
plotWells(2).bhp = true;
plotWells(5).bhp = true;
plotWells(12).qOs = true;
plotWells(12).qWs = true;

% To plot all injector bhp and all producer rates, uncomment the following
% lines:
if true
    for w = 1:numWells
        if w < 9
         plotWells(w).bhp = true;
        else
         plotWells(w).qOs = true;
         plotWells(w).qWs = true;
        end
    end
end
%%
%close all
plotSubIterations = false;
esmdaLegend = {'observations', 'truth', 'posterior mean', 'prior mean'};
if plotSubIterations
    esmdaLegend = {'observations', 'truth', 'posterior mean', 'ES-MDA it 3',...
                   'ES-MDA it 2', 'ES-MDA it 1', 'prior mean'}; %#ok
end

disp('Plotting in progress, please wait...');
ensemble.plotQoI('subplots', false, 'clearFigure', false, ...
    'cmapName', 'lines', ...
    'plotTruth', true, ...
    'subIterations', plotSubIterations, ...
    'observationIndices', observationIndices, ...
    'plotWells', plotWells, ... 
    'legend', esmdaLegend, ...
    'Position', [50 200 560 420], ...
    'savefig', true, ...
    'saveFolder', fullfile(historyMatchingDirectory, 'original_schedule'));
disp('Plotting completed');







%% Compare with original schedule
% Setup a new reference case using the original schedule, and then run the
% same schedule using the calibrated parameters

orginalReferenceExample = MRSTExample('egg_wo', 'realization', referenceEggRealization);
orginalReferenceExample.schedule = originalSchedule;

originalReferenceProblem = orginalReferenceExample.getPackedSimulationProblem('Directory', originalScheduleDirectory);

rerunOriginalSchedule = true;
if rerunOriginalSchedule
    clearPackedSimulatorOutput(originalReferenceProblem, 'prompt',  false);
end
simulatePackedProblem(originalReferenceProblem);

%% Create QoI using the original schedule as reference data
originalScheduleQoI = WellQoIHM('wellNames', wellNames, ...
                                'names', fieldNames);

originalScheduleQoI = originalScheduleQoI.validateQoI(originalReferenceProblem);
originalScheduleReferenceObservations = originalScheduleQoI.getQoI(originalReferenceProblem);

originalScheduleQoI.ResultHandler.dataPrefix = 'originalScheduleObservedQoI';

% Fill the the first resulthandler with the reference data
originalScheduleQoI.ResultHandler{1} = {originalScheduleReferenceObservations};

% Store the reference data in a dedicated result handler for comparing
% predictive performance of the calibrated models 
originalScheduleTruthResultHandler = ResultHandler('dataPrefix', 'trueQoI', ...
                                   'writeToDisk', originalScheduleQoI.ResultHandler.writeToDisk,...
                                   'dataDirectory', originalScheduleQoI.ResultHandler.dataDirectory, ...
                                   'dataFolder', originalScheduleQoI.ResultHandler.dataFolder, ...
                                   'cleardir', false);
originalScheduleTruthResultHandler{1} = {originalScheduleReferenceObservations};


%% Create ensemble using the calibrated samples
% Here, we also need to redefine the QoI to use the original schedule
% output as reference data so that the plotting becomes easy

calibratedSamples = ensemble.samples; %.getMeanSample();

% Quantity of interest
calibratedQoi = WellQoIHM('wellNames', wellNames, ...
                'names', fieldNames, ...
                'observationResultHandler', originalScheduleQoI.ResultHandler, ...
                'truthResultHandler', originalScheduleTruthResultHandler, ...
                'observationCov', obsStdDev.^2);

calibratedSimulationsDirectory = fullfile(topDirectory, ...
    ['original_schedule_upscale_ensemble_egg_', num2str(referenceEggRealization)]);

origSchedbaseUpscaledModel = MRSTExample('upscaled_coarse_network', ...
                                'partition', [6 6 1], ...
                                'referenceExample', orginalReferenceExample, ...
                                'plotCoarseModel', false);
%%
originalScheduleEnsemble = MRSTHistoryMatchingEnsemble(...
    origSchedbaseUpscaledModel, calibratedSamples, calibratedQoi, ...
    'directory', calibratedSimulationsDirectory, ...
    'simulationStrategy', 'spmd', ...
    'reset', true, ...
    'verboseSimulation', false);


%% Simulate 
originalScheduleEnsemble.simulateEnsembleMembers();
%originalScheduleEnsemble.updateHistoryMatchingInterval((1));


%% Plot the results
disp('Plotting in progress...');
originalScheduleEnsemble.plotQoI('subplots', false, 'clearFigure', false, ...
    'cmapName', 'lines', ...
    'plotTruth', true, ...
    'plotObservation', false, ...
    'observationIndices', [1], ...
    'plotWells', plotWells, ... 
    'legend', {'truth', ...
    'calibrated parameters'}, ...
    'Position', [700 200 560 420], ...
    'savefig', true, ...
    'saveFolder', fullfile(historyMatchingDirectory, 'original_schedule'));
    %'Position', [700 200 560 420]);
disp('Plotting completed');





%%
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
