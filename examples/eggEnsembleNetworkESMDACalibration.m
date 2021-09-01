%% Egg: Ensemble-based calibration of network models using ES-MDA
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
% Note that this example builds heavily on the ensemble module in MRST, and
% could just as well have been placed in that module. 


mrstModule add ad-core ad-blackoil mrst-gui ad-props ...
    example-suite incomp ensemble network-models diagnostics optimization

mrstVerbose off

% Some options for running the example
rerunReferenceModel = false; 
plotReferenceModel = true;
regenerateInitialEnsemble = false;

%% Define reference realization and where to store reference data
% Here, we choose which realization from the Egg ensemble we use as
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

referenceExample = MRSTExample('egg_wo', 'realization', referenceEggRealization);

% We consider only the first 48 time steps from the reference model, since
% most of the model dynamics play out during that periode.
numTotalTimesteps = 48;
%referenceExample.schedule = simpleSchedule(referenceExample.schedule.step.val(1:48), ...
%                                           'W', referenceExample.schedule.control.W);
referenceExample.schedule = simpleSchedule(referenceExample.schedule.step.val(5:53), ...
                                           'W', referenceExample.schedule.control.W);
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

wellNames = {referenceExample.schedule.control.W.name};
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
% The model is structured according to the MRSTExample class, which
% integrates directly into the ensemble module.
% For details, see the function injector_to_producer_network.m
baseNetworkModel = MRSTExample('injector_to_producer_network', ...
                               'cellsPerConnection', 10, ...
                               'referenceExample', referenceExample, ...
                               'plotNetwork', true);

%% Build initial network ensemble through flow diagnostics
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
% Since the flow diagnostics analalysis will only identify communications
% between a subset of all injector-producer pairs, the ensemble members
% will consist of networks of different sizes and shapes. We therefore map
% all the resulting networks back to the injectors-to-producers network.
% The closed connections are given parameters that ensure valid
% simulations, but that give negligible contribution to the fluid flows.
%
% For initial values for the well production indices, we use the mean
% contribution from the perforated cells of the full model as the value for
% the single perforated cell in the network model. An alternative is to use
% the sum of the perforated cells, and we compute (and store) both for
% convenience.

initializationType = 'fd_preprocessor';
postprocStateNumber = 20;

numInjectors = 8;
numProducers = 4;
numConnections = numInjectors*numProducers;
eggRealizations = [1:100];
ensembleSize = numel(eggRealizations);

% Minimal values representing "closed" connections
transMin = eps;
pvMin = 0.1;
    

initParamFolder = fullfile(mrstOutputDirectory(), ...
                         'eggEnsembleNetworkESMDACalibration_examples');
initParamFile = fullfile(initParamFolder, ...
                         [initializationType, '_parameters.mat']);

if exist(initParamFile, 'file') && ~regenerateInitialEnsemble
    % Read pre-computed initial ensemble
    preComputedEnsemble = load(initParamFile);
    transmissibilities = preComputedEnsemble.transmissibilities;
    poreVolumes = preComputedEnsemble.poreVolumes;
    wellProductionIndicesSum = preComputedEnsemble.wellProductionIndicesSum;
    wellProductionIndicesMean = preComputedEnsemble.wellProductionIndicesMean;
    
    
else
    % Use flow diagnostic analysis to obtain initial parameters
    
    % Start with values that correspond to "closed" connections only.
    transmissibilities = ones(numel(eggRealizations), numConnections).*transMin;
    poreVolumes = ones(numel(eggRealizations), numConnections).*pvMin;
    % The rows are organized as follows for the different connections:
    % 1 (inj) -  9 (prod)
    % 1 (inj) - 10 (prod)
    % 1 (inj) - 11 (prod)
    % 1 (inj) - 12 (prod)
    % 2 (inj) - 9 (prod)
    % 2 (inj) - 10 (prod)
    % ...
    % 8 (inj) - 11 (prod)
    % 8 (inj) - 12 (prod)
    
    wellProductionIndicesSum = zeros(numel(eggRealizations), numWells);
    wellProductionIndicesMean = zeros(numel(eggRealizations), numWells);
    
    for eggRealization = eggRealizations
        
        % Create the full Egg model realization
        ensembleExample = MRSTExample('egg_wo', 'realization', eggRealization);

        % We consider only the timesteps up until the step specified
        % in the variable 'postprocStateNumber'
        ensembleExample.schedule = simpleSchedule(ensembleExample.schedule.step.val(1:postprocStateNumber), ...
                                                   'W', ensembleExample.schedule.control.W);
       % Pack the problem 
       ensembleProblem = ensembleExample.getPackedSimulationProblem('Directory', ...
            fullfile(fullEnsembleDirectory, ['full_egg_', num2str(eggRealization)]));

        % Configure wells to only have one perforation and store sum of the
        % well production indices to use as prior.
        WtmpNetwork = ensembleExample.schedule.control.W;
        for i = 1:numel(WtmpNetwork)
            WtmpNetwork(i).cells = WtmpNetwork(i).cells(7);
            wellProductionIndicesSum(eggRealization, i) = sum(WtmpNetwork(i).WI);
            wellProductionIndicesMean(eggRealization, i) = mean(WtmpNetwork(i).WI);
        end
        
        % Network from flow diagnostic analysis
        switch initializationType
            case 'fd_preprocessor'
                % Flow diagnostic analysis directly on the geological model
                tmpNetwork = Network(WtmpNetwork, referenceExample.model.G, ...
                                     'type', initializationType,         ...
                                     'problem', ensembleProblem,        ...
                                     'flow_filter',1*stb/day);
            case 'fd_postprocessor'
                % Run the simulation 
                clearPackedSimulatorOutput(ensembleProblem, 'prompt',  false);
                simulatePackedProblem(ensembleProblem);
                
                % Then do postprocessing flow diagnostic analysis
                tmpNetwork = Network(WtmpNetwork, referenceExample.model.G, ...
                                     'type', initializationType,         ...
                                     'problem', referenceProblem,        ...
                                     'state_number',procStateNumber,     ...
                                     'flow_filter', 1*stb/day);
            otherwise
                error('\nNetwork of type %s is not implemented\n', initializationType);    
        end

        % Map the values from the smaller network to the full network
        for connectionID = 1:height(tmpNetwork.network.Edges)
            connection = tmpNetwork.network.Edges(connectionID,:);
            inj  = connection.EndNodes(1);
            prod = connection.EndNodes(2);
            %[inj, prod, (numProducers*(inj-1) + prod-numInjectors)]
            rowInFullNetwork = numProducers*(inj-1) + prod-numInjectors;
            
            transmissibilities(eggRealization, rowInFullNetwork) = connection.Transmissibility;
            poreVolumes(eggRealization, rowInFullNetwork) = connection.PoreVolume/baseNetworkModel.options.cellsPerConnection;
        end
        
        fprintf('Flow diagnostic analysis %d%% complete\n', eggRealization);
        
    end
    
    if ~exist(initParamFolder, 'dir')
        mkdir(initParamFolder);
    end
    save(initParamFile, 'transmissibilities', 'poreVolumes', ...
        'wellProductionIndicesSum', 'wellProductionIndicesMean');
end


%% Create sample object from initial ensemble
% We store the initial ensemble parameters in the appropriate Sample
% objects designed for history matching through the ensemble module. For 
% this, we use a class that are designed specifically for network models
% for holding the relevant parameters. We also specify maximum and minimum
% allowed values to avoid that simulations break down due to nonphysical
% configurations.

% The parameters transmissibility and pore volume belongs to the operators
% of a model.
operatorData = cell(ensembleSize, 1);
for i = 1:ensembleSize
    operatorData{i}.pv = poreVolumes(i, :);
    operatorData{i}.T = transmissibilities(i, :);
end
transMax = 1.5*max(max(transmissibilities));
pvMax = 1.5*max(max(poreVolumes));

connectionIndices = struct('faces', {baseNetworkModel.options.networkModel.Graph.Edges.Face_Indices}, ...
                           'cells', {baseNetworkModel.options.networkModel.Graph.Edges.Cell_Indices});
operatorSamples = NetworkOperatorSamplesHM('data', operatorData, ...
    'connectionIndices', connectionIndices, ...
    ... % parameter scaling 
    'pvScale', 1e4, 'TScale', 1e-9, ...
    'minPvValue', pvMin, 'maxPvValue', pvMax, ...
    'minTValue', transMin, 'maxTValue', transMax);

% Well production indices is considered a well parameter
wellSampleData = cell(ensembleSize, 1);

% We choose either the mean or the sum of the original well production
% indices
wellProductionIndices = wellProductionIndicesMean;

for i = 1:ensembleSize
    wellSampleData{i}.WI = wellProductionIndices(i, :);
end
minWI = 0.01*min(min(wellProductionIndices(:,:)));
maxWI = 7*max(max(wellProductionIndices(:,:)));

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
% choose parallel simulation strategy, etc.

ensemble = MRSTHistoryMatchingEnsemble(baseNetworkModel, samples, qoi, ...
    'alpha', [28/3 7 4 2], ...
    'directory', historyMatchingDirectory, ...
    'simulationStrategy', 'spmd', ...
    'maxWorkers', 8, ...
    'reset', true, ...
    'verbose', true, ...
    'verboseSimulation', false);

%% Run prior ensemble 
ensemble.simulateEnsembleMembers();

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
ensemble.simulateEnsembleMembers();

%% Plot prior and posterior ensemble results
close all
disp('Plotting in progress, please wait...');
ensemble.plotQoI('subplots', false, 'clearFigure', false, ...
    'cmapName', 'lines', ...
    'plotTruth', true, ...
    'subIterations', true, ...
    'observationIndices', observationIndices, ...
    'legend', {'observations', 'truth', 'posterior mean', 'ES-MDA it 3',...
               'ES-MDA it 2', 'ES-MDA it 1', 'prior mean'});
disp('Plotting completed');








