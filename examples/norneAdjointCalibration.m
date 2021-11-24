%% Norne: calibration of a GPSNet model
% This example demonstrates how to set up a GPSNet type reduced network
% model for a geostatistical realization of the Norne field model with
% simplified flow physics and well placement.
%
% To calibrate the model, we specify a training simulation with oscillatory
% well controls and use the Broyden–Fletcher–Goldfarb–Shanno (BFGS)
% algorithm. This is an iterative line-search method that gradually
% improves an approximation to the Hessian matrix of the mismatch function,
% obtained only from adjoint gradients via a generalized secant method.
%
% This example was first introduced in MRST 2021b.
mrstModule add ad-core ad-blackoil deckformat diagnostics...
               mrst-gui ad-props incomp optimization...
               network-models example-suite linearsolvers 

%% Setup 3D reference model
% The reference model is a single stochastic realization of the Norne field
% model taken from the example-suit module. Compared with the real field,
% this simulation case has simpler fluid description (an oil-water model)
% and an idealized field development plan consisting of a simple pattern of
% eleven vertical wells that run under constant bhp or rate controls.
%
% To train the GPSNet model, we set up a schedule that oscillates randomly
% around the prescribed liquid rates and pressure controls of the reference
% solution. This is done to generate a more general model that should
% provide accurate predictions in a sufficiently small region around the
% specified controls.

% True schedule, which we seek to reproduce
trueEx   = MRSTExample('norne_simple_wo');
trueCase = trueEx.getPackedSimulationProblem();
%clearPackedSimulatorOutput(trueProb)
simulatePackedProblem(trueCase);

[wellSolTrue, statesTrue] = getPackedSimulatorOutput(trueCase);
modelTrue     = trueEx.model;
scheduleTrue  = trueCase.SimulatorSetup.schedule;
WTrue         = scheduleTrue.control.W;

% Random schedule
trainEx   = makeRandomTraining(trueEx, false);
trainCase = trainEx.getPackedSimulationProblem();
%clearPackedSimulatorOutput(trainProb)
simulatePackedProblem(trainCase);

[wellSolTrain, statesTrain] = getPackedSimulatorOutput(trainCase);
modelTrain     = trainEx.model;
scheduleTrain  = trainCase.SimulatorSetup.schedule;
WTrain         = scheduleTrain.control.W;

% Plot
trueEx.plot(statesTrue,'step_index',numel(statesTrue))

plotWellSols({wellSolTrain, wellSolTrue}, ...
    {scheduleTrain.step.val, scheduleTrue.step.val},...
    'datasetnames',{'training','reference'}, ...
    'zoom', true, 'field', 'qWs', 'SelectedWells', 1:6);


%% Create the network
% We start by creating a network that connects injectors and producers.
% This network describes the possible 1D flow paths that each will have an
% associated pore volume and transmissibility.

% To create the network, we only need a single perforation from each well,
% which we somewhat arbitrarily pick to be the middle perforation.
Wnw = WTrue;
for i = 1:numel(Wnw)
    Wnw(i).cells = Wnw(i).cells(round(numel(Wnw(i).cells)/2));
end

% Select and create network
% The module offers different ways to specify the network topology. To
% choose different ones, you change the networkType variable
networkType = 'injectors_to_producers';
switch networkType
    case 'all_to_all'
        ntwrk =  Network(Wnw, modelTrue.G, 'type', networkType);
    case 'injectors_to_producers'
        ntwrk =  Network(Wnw, modelTrue.G, 'type', networkType, ...
                         'injectors', 1:6, 'producers', 7:11);
    case 'user_defined_edges'
        edges = [rldecode((1:6)',5*ones(6,1)), repmat((7:11)',6,1)];
        edges = sortrows([edges; edges],'ascend');
        ntwrk =  Network(Wnw, modelTrue.G, 'type', networkType, ...
                         'edges',edges);
    case 'fd_preprocessor'
        ntwrk =  Network(Wnw, modelTrue.G, 'type', networkType, ...
                         'problem', trueCase, 'flow_filter', 1*stb/day);
    case 'fd_postprocessor'
        ntwrk =  Network(Wnw, modelTrue.G, 'type', networkType, ...
                         'problem', trueCase,  ...
                         'state_number', 6,   ...
                         'flow_filter',1*stb/day);
    otherwise
        error('\nType of network: %s is not implemented\n', networkType);
end

%% Create the GPSNet
% Each edge in the network is subgridded and mapped onto a row in a
% rectangular Cartesian gird having the same number of rows as the number
% of network edges. The fluid model is copied from the reference model.
gravity off
gpsNet = NetworkModel(modelTrue, ntwrk, WTrain);

%% Plot the GPSNet model: network and simulation grid
fig1 = figure;
subplot(2,2,1)
G = trueEx.getVisualizationGrid();
plotCellData(G, trueEx.model.rock.poro,'EdgeColor','none');
plotWell(G, trueEx.schedule.control(1).W,'color','k','FontSize',10);
view(85,65); 
axis tight off; set(gca,'Clipping',false); zoom(1.2);

subplot(2,2,3), ax=gca;
cmap=ntwrk.plotNetwork();
axis tight off; set(ax,'Clipping',false); zoom(1.3); 

subplot(2,2,4), ax=gca;
ntwrk.plotNetwork('circle'); 
axis equal tight off; set(ax,'Clipping',false); zoom(1.3); 

subplot(2,2,2), ax=gca;
gpsNet.plotGrid([]);
colormap(ax,min(cmap+.2,1));
for i=1:5, ax.Children(i).Position(1:2) =  [4050,-10]; end
for i=6:11, ax.Children(i).Position(1:2) = [-125,-10]; end
set(ax.Children(1:11),'FontSize',8);
ax.Position([2 4]) = ax.Position([2 4])+[-.075 .1];

%% Specify training schedule and parameters to be matched
% We need to specify which of the parameters in the model we wish to match
% and any scaling and box limits that will applied to these in the BFGS
% optimization. These are specified through through the ModelParameter
% class from the optimization module. For the training we use the random
% schedule
trainSetup = gpsNetSimulationSetup(gpsNet, scheduleTrain);
predSetup  = gpsNetSimulationSetup(gpsNet, scheduleTrue);

% Parameter lumping
% In the GPSNet type of models, we only match a single pore volume and a
% single transmissibility for each network edge. We thus need a mapping
% between edges and all associated cells/faces in the Cartesian grid to
% describe the corresponding lumping of parameters in this grid.
[cellEdgeNo, cellIx] = gpsNet.getMapping('cells');
[faceEdgeNo, faceIx] = gpsNet.getMapping('faces');

config = {
    ...%name      include     scaling    boxlims     lumping     subset   relativeLimits
    'porevolume',       1,   'linear',       [],  cellEdgeNo,   cellIx,   [.001 4]
    'conntrans',        1,   'log',          [],          [],       [],   [.001 100]
    'transmissibility', 1,   'log'   ,       [],  faceEdgeNo,   faceIx,   [.001 100]};
prmsTrain = []; prmsTrue = [];
for k = 1:size(config,1)
    if config{k, 2} == 0, continue, end
    prmsTrain = addParameter(prmsTrain, trainSetup, ...
        'name',    config{k,1}, 'scaling', config{k,3}, ...
        'boxLims', config{k,4}, 'lumping', config{k,5}, ...
        'subset',  config{k,6}, 'relativeLimits',config{k,7});
    prmsTrue = addParameter(prmsTrue, predSetup, ...
        'name',    config{k,1}, 'scaling', config{k,3}, ...
        'boxLims', config{k,4}, 'lumping', config{k,5}, ...
        'subset',  config{k,6}, 'relativeLimits',config{k,7});
end


%% Define the mismatch function
% The mismatch function is defined as a function handle to a library
% function from the optimization module that computes the mismatch between
% a given simulation and a reference state. For an oil-water system, the
% match is computed based on three quantities (water/oil rate and bhp) and
% these must be given an associated weight.
weighting  = {'WaterRateWeight',  day/10000, ...
              'OilRateWeight',    day/20000, ...
              'BHPWeight',        1/(500*barsa)};
mismatchFn = @(model, states, schedule, states_ref, tt, tstep, state) ...
    matchObservedOW(model, states, schedule, states_ref,...
                   'computePartials', tt, 'tstep', tstep, weighting{:},...
                   'state',state,'from_states',false);


%% Evaluate the initial misfit in training data and prediction data
% To define the initial model, we extract the network parameters from the
% initial specification, scale them to the unit interval, and organize them
% into a vector that can be passed onto mismatchFn
pinit = gpsNet.getScaledParameterVector(trainSetup, prmsTrain, 0.5);
 
% Evaluate the initial mismatch for the training data
[misfitT0,~,wellSolsT0] = ...
    evaluateMatch(pinit,mismatchFn,trainSetup,prmsTrain,statesTrain,'Gradient','none');
[misfitE0,~,wellSolsE0] = ...
    evaluateMatch(pinit,mismatchFn,predSetup,prmsTrue,statesTrue,'Gradient','none');

%% Model calibration
% Calibrate the model using the BFGS method. This is a computationally
% expensive operation that may run for several hours if you choose a large
% number of iterations. Here, we therefore only apply 10 iterations. To get
% a good match, you should increase this number to 100+
objh = @(p)evaluateMatch(p,mismatchFn,trainSetup,prmsTrain,statesTrain);
[v, popt, history] = unitBoxBFGS(pinit, objh, 'objChangeTol', 1e-8, 'gradTol', 1e-5, ...
    'maxIt', 10, 'lbfgsStrategy', 'dynamic', 'lbfgsNum', 5, ...
    'outputHessian', true, 'logPlot', true);

%% Evaluate mismatch over the full simulation schedule 
[misfitE, ~,wellSolsE] = ...
    evaluateMatch(popt,mismatchFn,predSetup,prmsTrue,statesTrue,'Gradient','none');
[misfitT,~,wellSolsT] = ...
    evaluateMatch(popt,mismatchFn,trainSetup,prmsTrain,statesTrain,'Gradient','none');

%% Plot well responses for training data
trainSteps = scheduleTrain.step.val;
fh = plotWellSols({wellSolTrain,wellSolsT0,wellSolsT}, ...
    {trainSteps, trainSteps, trainSteps}, ...
    'datasetnames', {'train','init','match'}, 'zoom', true, ...
    'field', 'qWs', 'SelectedWells', 7:8);
set(fh, 'name','Norne: GPSNet training')


%% Plot well responses for prediction case
predSteps = scheduleTrue.step.val;
fh = plotWellSols({wellSolTrue,wellSolsE0, wellSolsE}, ...
    {predSteps, predSteps, predSteps}, ...
    'datasetnames', {'reference','initial','predicted'}, 'zoom', true, ...
    'field', 'qOs', 'SelectedWells', 7:8);
set(fh, 'name','Norne: GPSNet prediction')