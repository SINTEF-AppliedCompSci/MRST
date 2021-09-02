%% Norne: calibration of a GPSNet model
% This example demonstrates how to set up a GPSNet type reduced network
% model for a geostatistical realization of the Norne field model with
% simplified flow physics and well placement.
%
% To calibrate the model, we use the Broyden–Fletcher–Goldfarb–Shanno
% (BFGS) algorithm. This is an iterative line-search method that gradually
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
example = MRSTExample('norne_simple_wo');
problem = example.getPackedSimulationProblem();

% Simulate
simulatePackedProblem(problem);

[wellSolsRef, statesRef] = getPackedSimulatorOutput(problem);
modelRef    = example.model;
scheduleRef = problem.SimulatorSetup.schedule;
Wref        = scheduleRef.control.W;

% Plot
example.plot(statesRef,'step_index',numel(statesRef))

%% Create the network
% We start by creating a network that connects injectors and producers.
% This network describes the possible 1D flow paths that each will have an
% associated pore volume and transmissibility.

% To create the network, we only need a single perforation from each well,
% which we somewhat arbitrarily pick to be the middle perforation.
Wnw = Wref;
for i = 1:numel(Wnw)
    num_cells = numel(Wnw(i).cells);
    Wnw(i).cells = Wnw(i).cells(round(num_cells/2));
end

% Select and create network
% The module offers different ways to specify the network topology. To
% choose different ones, you change the networkType variable
networkType = 'fd_postprocessor';
switch networkType
    case 'all_to_all'
        ntwrk =  Network(Wnw, modelRef.G, 'type', networkType);
    case 'injectors_to_producers'
        ntwrk =  Network(Wnw, modelRef.G, 'type', networkType, ...
                         'injectors', 1:6, 'producers', 7:11);
    case 'fd_preprocessor'
        ntwrk =  Network(Wnw, modelRef.G, 'type', networkType, ...
                         'problem', problem, 'flow_filter', 1*stb/day);
    case 'fd_postprocessor'
        ntwrk =  Network(Wnw, modelRef.G, 'type', networkType, ...
                         'problem', problem,  ...
                         'state_number', 6,   ...
                         'flow_filter',1*stb/day);
    otherwise
        error('\nType of network: %s is not implemented\n', networkType);
end

% Plot the network.
% If based on flow diagnostics, we plot the network twice to show the
% relative magnitude of the associated transmissibilities and pore volumes.
figure; ntwrk.plotNetwork()

%% Create the data-driven model
% We subgrid each flow path into ten uniform cells and map the resulting
% network onto a rectangular Cartesian grid having the same number of rows
% as the number of flow paths.  The fluid model is copied from the model we
% seek to match.

% Grid and petrophysics
% The grid is set to have an aspect ratio of [5 1 1] and a volum that
% matches the bulk volume of the reservoir. The constant petrophysical
% properties are dummy values primarily used to compute an initial guess
% for matching parameters that have not be set by the network type.
L    = nthroot(sum(modelRef.operators.pv./modelRef.rock.poro)*25,3)  ;                            
G    = cartGrid([10, 1, numedges(ntwrk.network)], [L, L/5 ,L/5]*meter^3);
G    = computeGeometry(G);
rock = makeRock(G, 1000*milli*darcy, 0.2);

% Fluid model with endpoint scaling of relative permeabilities
gravity off
model = GenericBlackOilModel(G, rock, modelRef.fluid,'gas', false);
model = imposeRelpermScaling(model, 'SWL', .1, 'SWCR', .2, 'SWU', .9, ...
                             'SOWCR', .1, 'KRW', .9, 'KRO', .8);
model.OutputStateFunctions = {};

% Create the network model
% Here, we use a generic black-oil type model and it is thus important to
% validate the model to dynamically construct the correct number of phases
% and components and set up all the necessary state functions.
cellsPerPath = 10;
ntwrkModel = NetworkModel(model, cellsPerPath, ntwrk.network, Wref);
model      = ntwrkModel.model;
model      = model.validateModel();
W          = ntwrkModel.W;
state0     = initState(model.G, W, 100*barsa, [0,1]); 

% Simulation tolerances
model.toleranceCNV = 1e-6;
%model.FacilityModel.toleranceWellRate = 1e-4;
%model.FacilityModel.toleranceWellMS  = 1e-4;


%% Specify training schedule and parameters to be matched
% We need to specify which of the parameters in the model we wish to match
% and any scaling and box limits that will applied to these in the BFGS
% optimization. These are specified through through the ModelParameter
% class from the optimization module. For the training we use the first
% half of the reference simulation
trainSteps    = 1:round(numel(scheduleRef.step.val)/2);
trainSchedule = simpleSchedule(scheduleRef.step.val(trainSteps), 'W', W);
trainProb     = struct('model', model, 'schedule', trainSchedule, 'state0', state0);                               

% Parameter lumping
% In the GPSNet type of models, we only match a single pore volume and a
% single transmissibility for each network edge. We thus need a mapping
% between edges and cells/faces in the Cartesian grid to describe the
% corresponding lumping of parameters in this grid.
[cellEdgeNo, cellIx] = reorganizeIndices(ntwrkModel.Graph.Edges.Cell_Indices);
[faceEdgeNo, faceIx] = reorganizeIndices(ntwrkModel.Graph.Edges.Face_Indices);

nc =  model.G.cells.num;
nf =  numel(model.operators.T);

% select configuration for sensitivity computations
config = {
    ...%name      include     scaling    boxlims     lumping     subset   relativeLimits
    'porevolume',       1,   'linear',       [],  cellEdgeNo,   cellIx,   [.001 4]
    'conntrans',        1,   'log',          [],          [],       [],   [.001 100]
    'transmissibility', 1,   'log'   ,       [],  faceEdgeNo,   faceIx,   [.001 100]
    'swl',              1,   'linear',       [],  ones(nc,1),       [],   [.5 2]
    'swcr',             1,   'linear',       [],  ones(nc,1),       [],   [.5 2]
    'swu',              1,   'linear',       [],  ones(nc,1),       [],   [.5 2]
    'sowcr',            1,   'linear',       [],  ones(nc,1),       [],   [.5 2]
    'krw',              1,   'linear',       [],  ones(nc,1),       [],   [.5 2]
    'kro',              1,   'linear',       [],  ones(nc,1),       [],   [.5 2]
    'sw',               1,   'linear',   [0 .6],  cellEdgeNo,   cellIx,   [0 10]
    'pressure'          1,   'linear',       [],  cellEdgeNo,   cellIx,   [.1 4]};
parameters = [];
for k = 1:size(config,1)
    if config{k, 2} == 0, continue, end
    parameters = addParameter(parameters, trainProb, ...
        'name',    config{k,1}, 'scaling',       config{k,3}, ...
        'boxLims', config{k,4}, 'lumping',       config{k,5}, ...
        'subset',  config{k,6}, 'relativeLimits',config{k,7});
end


%% Define the mismatch function
% The mismatch function is defined as a function handle to a library
% function from the optimization module that computes the mismatch between
% a given simulation and a reference state. For an oil-water system, the
% match is computed based on three quantities (water/oil rate and bhp) and
% these must be given an associated weight.
weighting  = {'WaterRateWeight',  day/500, ...
              'OilRateWeight',    day/500, ...
              'BHPWeight',        1/(50*barsa)};
mismatchFn = @(model, states, schedule, states_ref, tt, tstep, state) ...
    matchObservedOW(model, states, schedule, states_ref,...
                   'computePartials', tt, 'tstep', tstep, weighting{:},...
                   'state',state,'from_states',false);


%% Set parameters defining the initial model
% To define the initial model, we extract the network parameters from the
% initial specification, scale them to the unit interval, and organize them
% into a vector that can be passed onto a library function from the
% optimization module that performs the forward simulation and evaluates
% the mismatch using mismatchFn.

% Extract and scale parameters
% Remember to overwrite pore volumes and transmissibilities if initial
% values for these have been computed by flow diagnostics
values = applyFunction(@(p)p.getParameterValue(trainProb), parameters);
values{2} =  0.5*values{2};
if any(strcmp(networkType,{'fd_preprocessor','fd_postprocessor'}))
    values{1} =  ntwrk.network.Edges.PoreVolume/cellsPerPath;
    values{3} =  ntwrk.network.Edges.Transmissibility;
end
for k = numel(values):-1:1    
    u{k} = parameters{k}.scale(values{k});
end
pvec0 = vertcat(u{:});
 
%{
% Perform a forward simulation and evaluate the mismatch
[misfitVal0,~,wellSols0,states0] = ...
    evaluateMatch(pvec0,mismatchFn,trainProb,parameters,statesRef,'Gradient','none');          

% Plot well curves
fh = plotWellSols({wellSolsRef,wellSols0}, ...
    {scheduleRef.step.val,trainProb.schedule.step.val}, ...
    'datasetnames', {'reference','initial'}, 'zoom', true, ...
    'field', 'qOs', 'SelectedWells', 7);
set(fh, 'name','Norne')
legend('reference model','initial DD model')
%}

%% Model calibration
% Calibrate the model using the BFGS method
objh = @(p)evaluateMatch(p,mismatchFn,trainProb,parameters,statesRef);

[v, pvecOpt, history] = unitBoxBFGS(pvec0, objh,'objChangeTol',  1e-8, ...
    'maxIt', 30, 'lbfgsStrategy', 'dynamic', 'lbfgsNum', 5);

%% Evaluate mismatch over the full simulation schedule
prob = trainProb;
prob.schedule = simpleSchedule(scheduleRef.step.val, 'W', W);
 
[~,~,wellSolsOpt] = evaluateMatch(pvecOpt,mismatchFn,prob,parameters, statesRef,'Gradient','none');
[~,~,wellSols0]   = evaluateMatch(pvec0,mismatchFn,prob,parameters, statesRef,'Gradient','none');

%% Plot well curves
fh = plotWellSols({wellSolsRef,wellSols0,wellSolsOpt}, ...
    {scheduleRef.step.val,prob.schedule.step.val,prob.schedule.step.val}, ...
    'datasetnames', {'reference','initial','optimized'}, 'zoom', true, ...
    'field', 'qOs', 'SelectedWells', 7);
set(fh, 'name','Norne')
legend('reference model','initial DD model','optimize DD model')
                                       