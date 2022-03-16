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
    test-suite incomp ensemble optimization upscaling coarsegrid

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

referenceCase = TestCase('egg_wo', 'realization', referenceEggRealization);

% We consider only the first 48 time steps from the reference model, since
% most of the model dynamics play out during that periode.
numTotalTimesteps = 48;

referenceCase.schedule = simpleSchedule(referenceCase.schedule.step.val(1:48), ...
                                           'W', referenceCase.schedule.control.W);
    
%% Create random schedule

originalSchedule = referenceCase.schedule;
schedule = referenceCase.schedule;

% Define new schedule after every fourth timestep
schedule.step.control= ceil(0.25*(1:numel(schedule.step.val))');

rng(0)
referenceCase.schedule = ensembleModulePerturbSchedule(schedule);

%% Pack and run reference problem
referenceProblem = referenceCase.getPackedSimulationProblem('Directory', referenceDirectory);

if rerunReferenceModel
    clearPackedSimulatorOutput(referenceProblem, 'prompt',  false);
end
simulatePackedProblem(referenceProblem);

if plotReferenceModel
    [refWellSols, refStates, refReports] = getPackedSimulatorOutput(referenceProblem);
    referenceCase.plot(refStates);
    plotWellSols(refWellSols);
end


%% Define the coarse upscaled model
% We create a rich network in the form of a coarse model which we will 
% calibrate according to the reference data. 
%
% The model is structured according to the MRSTExample class, which
% integrates directly into the esemble module.
% For details, see the function upscaled_coarse_network

baseUpscaledModel = TestCase('upscaled_coarse_network', ...
                             'partition', [6 6 1], ...
                             'referenceCase', referenceCase, ...
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
%
% We store the initial ensemble parameters in the appropriate Sample
% objects designed for history matching through the ensemble module. Both
% pore volumes and transmissibilities belong to the operator sample class,
% whereas well productivity index is found in well samples. We also specify 
% maximum and minimum allowed values to avoid that simulations break down 
% due to nonphysical configurations.
%
% All this is done in the following utility function:

eggRealizations = [1:10];
%eggRealizations = [1:100];

if (numel(eggRealizations) < 50)
    warning('A small ensemble can be used to investigate what is going on in the example, but you should use closer to 100 ensemble members for good results.');
end


samples = priorSamplesForUpscaledEggModel(baseUpscaledModel, ...
                                          referenceCase, ...
                                          regenerateInitialEnsemble, ...
                                          'eggRealizations', eggRealizations);



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
qoi = WellQoIHM('wellNames', {referenceCase.schedule.control(1).W.name}, ...
                'names', {'qOs', 'qWs', 'bhp'}, ...
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
    'observationProblem', referenceProblem, ...
    'perturbObservations', false, ...
    'alpha', [28/3 7 4 2], ...
    'directory', historyMatchingDirectory, ...
    'simulationStrategy', 'background', ...
    'maxWorkers', 4, ...
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


%% Do ensemble-based calibration using history matching algorithms
% This function performs ES-MDA updates with the relaxation parameters
% given as 'alpha' to the ensemble, running intermediate ensemble 
% simulations as it goes.
ensemble.doHistoryMatching()

%% Run the posterior ensemble using the calibrated parameters
% We specify that we want to run the posterior for the complete simulation
% time span.
ensemble.simulateEnsembleMembers('progressTitle', 'Simulating posterior ensemble');


%% Plot prior and posterior ensemble results
%
% Create a boolean array of structs that reflects the QoI object, where
% each boolean value reflects which value from which well will be plotted.
% Here, we chose to plot the bhp from the injectors and the production
% rates from the producers.
numWells = numel(ensemble.qoi.wellNames);
injectors = 1:8;
producers = 9:12;

plotWells = repmat( struct('qOs', false, 'qWs', false, 'bhp', false), 1, 12);
[plotWells(producers).qOs] = deal(true);
[plotWells(producers).qWs] = deal(true);
[plotWells(injectors).bhp] = deal(true);

% Set this variable to true to plot intermediate ES-MDA iteration results.
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
    'plotWells', plotWells, ... 
    'legend', esmdaLegend, ...
    'Position', [50 200 560 420], ...
    'savefig', true, ...
    'saveFolder', fullfile(historyMatchingDirectory, 'original_schedule'));
disp('Plotting completed');



%% Compare with original schedule
% Setup a new reference case using the original schedule, and then run the
% same schedule using the calibrated parameters

orginalReferenceCase = TestCase('egg_wo', 'realization', referenceEggRealization);
orginalReferenceCase.schedule = originalSchedule;

originalReferenceProblem = orginalReferenceCase.getPackedSimulationProblem('Directory', originalScheduleDirectory);

rerunOriginalSchedule = true;
if rerunOriginalSchedule
    clearPackedSimulatorOutput(originalReferenceProblem, 'prompt',  false);
end
simulatePackedProblem(originalReferenceProblem);



%% Create ensemble using the calibrated samples
% Here, we use the calibrated samples from the previous ensemble.
% We also rebuild a base model using the original schedule.

calibratedSimulationsDirectory = fullfile(topDirectory, ...
    ['original_schedule_upscale_ensemble_egg_', num2str(referenceEggRealization)]);

origSchedbaseUpscaledModel = TestCase('upscaled_coarse_network', ...
                                'partition', [6 6 1], ...
                                'referenceCase', orginalReferenceCase, ...
                                'plotCoarseModel', false);

originalScheduleEnsemble = MRSTHistoryMatchingEnsemble(...
    origSchedbaseUpscaledModel, ensemble.samples, qoi, ...
    'observationProblem', originalReferenceProblem, ...
    'perturbObservations', false, ...
    'directory', calibratedSimulationsDirectory, ...
    'simulationStrategy', 'spmd', ...
    'reset', true, ...
    'verboseSimulation', false);


%% Simulate upscaled 
% Using the original schedule and the calibrated parameters
originalScheduleEnsemble.simulateEnsembleMembers();


%% Plot the results
disp('Plotting in progress...');
originalScheduleEnsemble.plotQoI('subplots', false, 'clearFigure', false, ...
    'cmapName', 'lines', ...
    'plotTruth', true, ...
    'plotObservation', false, ...
    'plotWells', plotWells, ... 
    'legend', {'truth', 'calibrated parameters'}, ...
    'Position', [700 200 560 420], ...
    'savefig', true, ...
    'saveFolder', fullfile(historyMatchingDirectory, 'original_schedule'));
disp('Plotting completed');





%%
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
