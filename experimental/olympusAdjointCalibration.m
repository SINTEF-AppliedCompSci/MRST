%% Norne: calibration of a GPSNet model
% This example demonstrates how to set up a GPSNet type reduced network
% model for a geostatistical realization of the Norne field model with
% simplified flow physics and well placement.
%
% To calibrate the model, we use the Broyden–Fletcher–Goldfarb–Shanno
% (BFGS) algorithm. This is an iterative line-search method that gradually
% improves an approximation to the Hessian matrix of the mismatch function,
% obtained only from adjoint gradients via a generalized secant method.
mrstModule add ad-core ad-blackoil deckformat diagnostics...
               mrst-gui ad-props incomp optimization...
               network-models example-suite linearsolvers 
         
%% Setup 3D reference model
% The reference model is a single stochastic realization of the Olympus
% field model (the first ensemble member). Because MRST does not support
% the TUNING keyword in ECLIPSE input decks, we set the schedule manually.
example = MRSTExample('olympus_field_wo');
%example.plot(example.model.rock);

% Set up packed simulation
problem = example.getPackedSimulationProblem();
%clearPackedSimulatorOutput(problem)

% Fix schedule to use more reasonable time steps
% The input conversion only understands the DATES keyword and hence ends up
% with one-year time steps. We split each of these time steps in twelve to
% get closer to the maximum time step of 30 days specified by the TUNING
% keyword.
schedule = example.schedule;
ix  = repmat(1:numel(schedule.step.val), [12 1]);
schedule.step.val = schedule.step.val(ix(:))/12;
schedule.step.control = schedule.step.control(ix(:));
problem.SimulatorSetup.schedule = schedule;
example.schedule = schedule;

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
% associated pore volume and transmissibility. Because the original model
% has multiple fluid regions, we cannot use flow diagnostics preprocessing,
% since the incompressible flow solvers this relies on do not yet handle
% such fluid object. 
Wnw = Wref;
networkType = 'fd_postprocessor';
%networkType = 'injectors_to_producers';
switch networkType
    case 'all_to_all'
        for i = 1:numel(Wnw)
            num_cells = numel(Wnw(i).cells);
            Wnw(i).cells = Wnw(i).cells(round(num_cells/2));
        end
        ntwrk =  Network(Wnw, modelRef.G, 'type', networkType);
    case 'injectors_to_producers'
        for i = 1:numel(Wnw)
            num_cells = numel(Wnw(i).cells);
            %Wnw(i).cells = Wnw(i).cells([1 round(num_cells/2) num_cells]);
            Wnw(i).cells = Wnw(i).cells(round(num_cells/2));
        end
        ntwrk =  Network(Wnw, modelRef.G, 'type', networkType, ...
                         'injectors', 1:7, 'producers', 8:18);
    case 'fd_postprocessor'
        for i = 1:numel(Wnw)
            num_cells = numel(Wnw(i).cells);
            Wnw(i).cells = Wnw(i).cells(round(num_cells/2));
        end
        ntwrk =  Network(Wnw, modelRef.G, 'type', networkType, ...
                         'problem', problem,                   ...
                         'flow_filter', 1*stb/day,             ...
                         'state_number', numel(statesRef));
    otherwise
        error('\nNetwork of type %s is not implemented\n', networkType);           
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
L    = nthroot(sum(modelRef.operators.pv./modelRef.rock.poro)*25,3);
G    = cartGrid([10, 1, numedges(ntwrk.network)], [L, L/5 ,L/5]*meter^3);
G    = computeGeometry(G);
rock = makeRock(G,500*milli*darcy, 0.2);

% We replace the fluid model, which has different fluid regions, by a a
% single fluid having quadratic relative permeabilities. The petrophysical
% behavior of the resulting black-oil model can then only be adjusted
% through its endpoint scaling.
fluid = initSimpleADIFluid('phases', 'WO',... 
                           'mu' , [.3, 3]*centi*poise,...
                           'rho', [1014, 859]*kilogram/meter^3, ...
                           'n', [2 2]);
fluid.krPts  = struct('w', [0 0 1 1], 'ow', [0 0 1 1]);
scaling      = {'SWL', .1, 'SWCR', .2, 'SWU', .9, ...
                'SOWCR', .1, 'KRW', .9, 'KRO', .8};

gravity off
model = GenericBlackOilModel(G, rock, fluid,'gas', false);
model = imposeRelpermScaling(model, scaling{:});
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
state0     = initState(G, W , 215*barsa,[0, 1]); 

% Simulation tolerances
model.toleranceCNV = 1e-6;

%% Specify training schedule and parameters to be matched
% We need to specify which of the parameters in the model we wish to match
% and any scaling and box limits that will applied to these in the BFGS
% optimization. These are specified through through the ModelParameter
% class from the optimization module. For the training we use the first
% half of the reference simulation
trainSteps = 1:185;
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
    'porevolume',       1,   'linear',       [],  cellEdgeNo,   cellIx,   [.001 2]
    'conntrans',        1,   'log',          [],          [],       [],   [.001 10]
    'transmissibility', 1,   'log'   ,       [],  faceEdgeNo,   faceIx,   [.001 4]
    'swl',              1,   'linear',       [],  ones(nc,1),       [],   [.5 2]
    'swcr',             1,   'linear',       [],  ones(nc,1),       [],   [.5 2]
    'swu',              1,   'linear',       [],  ones(nc,1),       [],   [.5 2]
    'sowcr',            1,   'linear',       [],  ones(nc,1),       [],   [.5 2]
    'krw',              1,   'linear',       [],  ones(nc,1),       [],   [.5 2]
    'kro',              1,   'linear',       [],  ones(nc,1),       [],   [.5 2]
    'sw',               1,   'linear',  [0 .99],  cellEdgeNo,   cellIx,   [0 9]
    'pressure'          1,   'linear',       [],  cellEdgeNo,   cellIx,   [.001 10]};
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
% these must be given an associated weight. Because the wells are
% controlled by bhp, we only match rates here.
weighting =  {'WaterRateWeight',  day/100, ...
              'OilRateWeight',    day/10};    
mismatchFn = @(model, states, schedule, states_ref, tt, tstep, state) ...
    matchObservedOW(model, states, schedule, states_ref,...
                    'computePartials', tt, 'tstep', tstep, weighting{:},...
                    'state', state, 'from_states', false);
          
%% Set parameters defining the initial model
% To define the initial model, we extract the network parameters from the
% initial specification, scale them to the unit interval, and organize them
% into a vector that can be passed onto a library function from the
% optimization module that performs the forward simulation and evaluates
% the mismatch using mismatchFn.
values = applyFunction(@(p)p.getParameterValue(trainProb), parameters);
% scale values
if any(strcmp(networkType,{'fd_preprocessor','fd_postprocessor'}))
     values{1} =  pv/10;
     values{3} =  TT;
else
  rng(12345)
  values{10} =  rand(size(values{10}));  
end
for k = numel(values):-1:1    
    u{k} = parameters{k}.scale(values{k});
end
pvec0 = vertcat(u{:});  
 
%% Model calibration with BFGS
objh = @(p)evaluateMatch(p,mismatchFn,trainProb,parameters,statesRef);

[v, pvecOpt, history] = ...
    unitBoxBFGS(pvec0, objh,'objChangeTol',  1e-8, 'maxIt', 2, ...
                'lbfgsStrategy', 'dynamic', 'lbfgsNum', 5);

%% Simulating base case and calibrated model for comparison
prob = trainProb;
prob.schedule = simpleSchedule(scheduleRef.step.val, 'W', W);
 
[~,~,wellSolsOpt] = evaluateMatch(pvecOpt,mismatchFn,prob,parameters,statesRef,'Gradient','none');
[~,~,wellSols0]   = evaluateMatch(pvec0,  mismatchFn,prob,parameters,statesRef,'Gradient','none');


%% Plot well curves
fh = plotWellSols({wellSolsRef,wellSols0,wellSolsOpt}, ...
    {scheduleRef.step.val,prob.schedule.step.val,prob.schedule.step.val}, ...
    'datasetnames', {'reference','initial','optimized'}, 'zoom', true, ...
    'field', 'qOs', 'SelectedWells', 8);
set(fh, 'name','Olympus network model')
legend('reference','initial network model','calibrated network model')
