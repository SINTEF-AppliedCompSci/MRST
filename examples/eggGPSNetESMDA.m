%% Ensemble-based calibration of network models of Egg using ES-MDA
% The aim of this example is to calibrate an ensemble of GPSNet type
% network models to match reference well data from a full Egg model. To do
% so, we generate network models using flow diagnostics on each of the
% realizatons from the Egg ensemble. This ensemble will then act as our
% prior, and we use ES-MDA to obtain a posterior parameter distribution. We
% consider the mean of these parameter as the calibrated model parameters.
% 
% The model is a water-oil reservoir, and we calibrate pore volumes,
% transmissibilities, and well production indices.
% 
% Note that this example builds heavily on the ensemble module in MRST.
%
% This example was first introduced in MRST 2021b.


mrstModule add ad-core ad-blackoil mrst-gui ad-props ...
    test-suite incomp ensemble network-models diagnostics optimization

mrstVerbose off

warning('This example requires some time to run.')

% Some options for running the example
rerunReferenceModel = false; 
plotReferenceModel = true;
regenerateInitialEnsemble = false;
forceRun = true; 


%% Define reference realization and where to store reference data
% We choose which realization from the Egg ensemble we use as
% reference data, and define where to store the data generated from this
% example.

referenceEggRealization = 0;

topDirectory = fullfile(mrstOutputDirectory(), 'network_models', 'eggEnsembleESMDA');
referenceDirectory = fullfile(topDirectory, 'truth', ...
    ['egg_realization_', num2str(referenceEggRealization)]);
fullEnsembleDirectory = fullfile(topDirectory, 'full_models');
historyMatchingDirectory = fullfile(topDirectory, ...
    ['network_ensemble_egg_', num2str(referenceEggRealization)]);

%% Setup full 3D reference model
% The reference model is one of the realizations from the Egg ensemble. It
% consists of a 3D oil-water model of a highly chanalized reservoir, driven
% by four producer wells and eight injectors.

referenceCase = TestCase('egg_wo', 'realization', referenceEggRealization);

% Most of the dynamics in the Egg model plays out during the first half or 
% so of the total time span. We therefore consider just a subset of the 
% schedule from the reference model.
referenceCase.schedule = simpleSchedule(referenceCase.schedule.step.val(5:53), ...
                                           'W', referenceCase.schedule.control.W);
numTotalTimesteps = numel(referenceCase.schedule.step.val);

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


%% Define the network 
% We create a network model which we will calibrate according to the
% reference data. Since we at this point do not have any information about
% the flow paths in the model, we define our network model to connect all
% injectors with all producers with single 1D flow paths. 
%
% We split each flow path into ten cells and define a computational 2D grid
% in which each flow path is represented with a single row. We use the same
% fluid model as in the full reference model.
%
% Since the pore volume, transmissibilities, and well production indices
% will be calibrated through an ensemble, the corresponding values in the
% reference model defined here will not be used. We therefore assign dummy
% values for these parameters.
%
% The model is structured according to the TestCase class, which
% integrates directly into the ensemble module.
% For details, see the function injector_to_producer_network.m

baseNetworkModel = TestCase('injector_to_producer_network', ...
                            'cellsPerConnection', 10, ...
                            'referenceCase', referenceCase, ...
                            'plotNetwork', true);

%% Create sample object with initial networks obtained through flow diagnostics
% We obtain the prior distribution of the model parameters from building
% network models of realizations 1 to 100 from the full models in the Egg
% ensemble. With eight injectors and four producers we have a network
% consisting of 32 connections, and we assign homogeneous pore volumes and
% transmissibilities for each connection. The network models will therefore
% have 32 parameters of each of the two parameters.
%
% Here, we can choose between using flow diagnostic analysis of
% either the 3D geological model (preproccessing) or the fine-scale
% reference simulation (postprocessing).
%
% Since the flow diagnostics analysis will only identify communications
% between a subset of all injector-producer pairs, the ensemble members
% will consist of networks of different sizes and shapes. We therefore map
% all the resulting networks back to the injectors-to-producers network.
% The closed connections are given parameters that ensure valid
% simulations, but that give negligible contribution to the fluid flow.
%
% For initial values for the well production indices, we use the mean
% contribution from the perforated cells of the full model as the value for
% the single perforated cell in the network model. An alternative is to use
% the sum of the perforated cells, and we compute (and store) both for
% convenience.
%
% We store the initial ensemble parameters in the appropriate Sample
% objects designed for history matching from the ensemble module. For 
% this, we use a class that are designed specifically for network models
% for holding the relevant parameters. We also specify maximum and minimum
% allowed values to avoid that simulations break down due to nonphysical
% configurations.
%
% For details, see createEggSamplesFromFlowDiagnostics.m

initializationType = 'fd_preprocessor';
%initializationType = 'fd_postprocessor';
eggRealizations = [1:10];
%eggRealizations = [1:100];

if (numel(eggRealizations) < 50)
    warning('A small ensemble can be used to investigate what is going on in the example, but you should use closer to 100 ensemble members for good results.');
end

samples = createEggSamplesFromFlowDiagnostics(...
    baseNetworkModel, eggRealizations, initializationType, ...
    'regenerateInitialEnsemble', regenerateInitialEnsemble, ...
    'fullEnsembleDirectory', fullEnsembleDirectory, ...
    'WItype', 'mean' ...
);

%% Create QoI object and define observation uncertainty
% The quantity of interest for an ensemble simulation defines which data we
% will keep after simulating each ensemble member.
% Here, this will be the production rates and bhp.
% This will also correspond to the data that we use in the history matching
% for calibration.
%
% In history matching, the observations are always assumed to contain
% some noise or uncertainty, and history matching algorithms, such as
% ES-MDA, is constructed to calibrate model parameters with respect to the
% uncertainty of the observations. In this example, we consider the
% reference data to be exact, but we still need to give a measure of the
% observation uncertainty for the algorithm to make sense. 

% Observation uncertainty
obsStdDevFlux = 50*stb()/day();
obsStdDevBhp  = 1*barsa();
obsStdDev = [obsStdDevFlux, obsStdDevFlux, obsStdDevBhp];

wellNames = {referenceCase.schedule.control.W.name};

% Quantity of interest
qoi = WellQoIHM('wellNames', wellNames, ...
                'names', {'qOs', 'qWs', 'bhp'}, ...
                'observationCov', obsStdDev.^2);
            
%% Create the ensemble object
% Here, we gather all the information we have made until now. The ensemble
% consists of the baseExample that defines everything the ensemble members
% have in common (geological model, grid, fluid, simulation methods, and so
% on). We also specify the alpha parameters for ES-MDA, choose parallel
% simulation strategy, etc. More information about input parameters can be 
% found in the base class MRSTEnsemble.m.

% The maxWorkers argument determines the maximum number of simultaneous 
% processes which will be used to simulate the ensemble. Increasing this
% number will increase the number of processes used, but care must be taken
% to ensure that your computer is capable of handling this number of 
% processes. Specifying too high a number will result in high memory usage 
% which could drastically slow down your computer.

ensemble = MRSTHistoryMatchingEnsemble(baseNetworkModel, samples, qoi, ...
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

% Running simulations for all ensemble members will take a long time,
% prompt the user before beginning simulations. 

ensemble.simulateEnsembleMembers('progressTitle', 'Simulating prior ensemble');

%% Configure the schedule used during history matching
% During the calibration, we only use every second observation, and we only
% consider the first half of the simulation time span. The second half will
% be used to look at predictive qualities of the results.

totalNumberOfTimesteps = numel(ensemble.originalSchedule.step.val);
observationIndices = (2:2:floor(totalNumberOfTimesteps/2));
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
numWells = numel(ensemble.setup.schedule.control.W);
injectors = find([ensemble.setup.schedule.control.W.sign] ==   1);
producers = find([ensemble.setup.schedule.control.W.sign] ==  -1);

plotWells = repmat( struct('qOs', false, 'qWs', false, 'bhp', false), 1, numWells);

% Plot bhp from wells 2 and 5, and the production rates from well 12 (PROD4)
plotWells(2).bhp = true;
plotWells(5).bhp = true;
plotWells(12).qOs = true;
plotWells(12).qWs = true;

% Plot bhp of all injectors and production rates from all producers
% [plotWells(producers).qOs] = deal(true);
% [plotWells(producers).qWs] = deal(true);
% [plotWells(injectors).bhp] = deal(true);


close all
disp('Plotting in progress, please wait...');
ensemble.plotQoI('subplots', false, 'clearFigure', false, ...
    'cmapName', 'lines', ...
    'plotTruth', true, ...
    'subIterations', true, ...
    'observationIndices', observationIndices, ...
    'plotWells', plotWells, ... 
    'legend', {'observations', 'truth', 'posterior mean', 'ES-MDA it 3',...
               'ES-MDA it 2', 'ES-MDA it 1', 'prior mean'});
disp('Plotting completed');







%% Running model using the mean parameter values
% (optional)
%  Often with model calibration, we wish to have a model with a single set
%  of parameters at the end. We therefore compute the mean of the prior and
%  the posterior and run the network model to compare the usefulness of the
%  mean parameter distributions.

meanPriorSample = samples.getMeanSample();
meanPosteriorSample = ensemble.samples.getMeanSample();

%% Define new ensemble object 
%  Even though we now are interested in single simulations of the prior and
%  posterior means, respectively, we use the history matching ensemble
%  class to easily run the simulation and get the plotting.
%  We also reuse the QoI object from before.

meanSimulationsDirectory = fullfile(topDirectory, ...
    ['means_network_ensemble_egg_', num2str(referenceEggRealization)]);

meanEnsemble = MRSTHistoryMatchingEnsemble(baseNetworkModel, meanPriorSample, qoi, ...
    'observationProblem', referenceProblem, ...
    'perturbObservations', false, ...
    'directory', meanSimulationsDirectory, ...
    'simulationStrategy', 'serial', ...
    'reset', true, ...
    'verboseSimulation', true);

% Run simulation with the mean of the prior parameters
meanEnsemble.simulateEnsembleMembers();

%% Manually update the ensemble samples to the posterior mean
%  When doing history matching, the ensemble class computes a new set of
%  parameters and updates the ensemble samples behind the scenes. Here, we
%  already have the new sample and update the ensemble samples directly.
meanEnsemble.updateHistoryMatchingIterations(meanPosteriorSample);

% Run simulation with the mean of the posterior parameters
meanEnsemble.simulateEnsembleMembers();

%% Plot the results using the prior and posterior mean parameters
meanEnsemble.plotQoI('subplots', false, 'clearFigure', false, ...
    'cmapName', 'lines', ...
    'plotTruth', true, ...
    'observationIndices', observationIndices, ...
    'plotWells', plotWells, ... 
    'legend', {'observations', 'truth', ...
    'posterior parameters mean', 'prior parameters mean'}, ...
    'Position', [700 200 560 420]);
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
