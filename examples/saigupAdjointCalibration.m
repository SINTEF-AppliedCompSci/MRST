%% Calibration of a GPSNet model for SAIGUP using adjoints
% This example demonstrates how to set up a GPSNet type reduced network
% model based on injector-productor connections for a reservoir with
% realistic geology, but somewhat simplified flow physics and well
% placement. To this end, we use a realization from the SAIGUP project,
% which describes a shallow marine reservoir typical of the North Sea using
% a 40x120x20 corner-point grid with 78200 active cells.
%
% To calibrate the model, we use the Broyden–Fletcher–Goldfarb–Shanno
% (BFGS) algorithm. This is an iterative line-search method that gradually
% improves an approximation to the Hessian matrix of the mismatch function,
% obtained only from adjoint gradients via a generalized secant method.
mrstModule add ad-core ad-blackoil deckformat diagnostics mrst-gui ...
    ad-props incomp optimization example-suite linearsolvers 

%% Setup 3D reference model
% Our reference model is the standard SAIGUP test case from the
% example-suite tutorial. This specifies an oil-water fluid model with a
% 5:1 viscosity ratio, quadratic relative permeabilities, and slight fluid
% compressibilities. Oil and water are initially separated by a sharp
% interface and the reservoir is produced using six vertical producers
% supported by eight vertical injectors around the perimeter of the model.
% If simulation results are not already stored for the model, we perform a
% full 3D simulation first.

example = MRSTExample('saigup_wo');
% example.plot(example.model.rock);

problem = example.getPackedSimulationProblem();
simulatePackedProblem(problem);

[wellSolsRef, statesRef] = getPackedSimulatorOutput(problem);
modelRef    = example.model;
scheduleRef = problem.SimulatorSetup.schedule;
Wref        = scheduleRef.control.W;
% example.plot(states);

%% Create the network
% We start by creating a network that connects injectors and producers.
% This network describes the possible 1D flow paths that each will have an
% associated pore volume and transmissibility.

% To create the network, we only need a single perforation from each well,
% which we somewhat arbitrarily pick to be the 7th from the top.
Wnw = Wref;
for i = 1:numel(Wnw)
    Wnw(i).cells = Wnw(i).cells(7);
end

% Select type and create network 
% The module supplies different ways to do this: all to all, injectors to
% producers (and vice versa), connections specified manually by the user,
% connections determined by a flow diagnostics analysis of the 3D
% geological model, or connections based on flow diagnostics analysis of
% the fine-scale reference simulation.

%networkType = 'all_to_all';
networkType  = 'injectors_to_producers';
%networkType = 'fd_preprocessor';
%networkType = 'fd_postprocessor';
switch networkType
    case 'all_to_all'
        ntwrk =  Network(Wnw, modelRef.G, ...
                         'type', networkType);
    case 'injectors_to_producers'
        ntwrk =  Network(Wnw, modelRef.G,      ...
                         'type', networkType,  ...
                         'injectors', 1:8,     ...
                         'producers', 9:17);
    case 'fd_preprocessor'
        ntwrk =  Network(Wnw,modelRef.G,       ...
                         'type', networkType,  ...
                         'problem', problem,   ...
                         'flow_filter',1*stb/day);
    case 'fd_postprocessor'
        ntwrk =  Network(Wnw,modelRef.G,       ...
                         'type', networkType,  ...
                         'problem', problem,   ...
                         state_number',20,     ...
                         'flow_filter', 1*stb/day);
    otherwise
        error('\nType of network: %s is not implemented\n', networkType);             
end
                     
% Plot the network.
% If based on flow diagnostics, we plot the network twice to show the
% relative magnitude of the associated transmissibilities and pore volumes.
if any(strcmp(networkType,{'fd_preprocessor','fd_postprocessor'}))
    TT = ntwrk.network.Edges.Transmissibility;
    pv = ntwrk.network.Edges.PoreVolume;

    figure, subplot(1,2,1);
    ntwrk.plotNetwork('NetworkLineWidth',10*TT/max(TT));
    title('Transmissibility');
    axis off

    subplot(1,2,2);
    ntwrk.plotNetwork('NetworkLineWidth',10*pv/max(pv));
    title('PoreVolume');
    axis off;
else
    figure
    ntwrk.plotNetwork()
    view(-80,80), axis off
end

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
rock = makeRock(G, 200*milli*darcy, 0.12);

% Fluid model with endpoint scaling of relative permeabilities
gravity off
fluid = modelRef.fluid;
model = GenericBlackOilModel(G, rock, fluid,'gas', false);
model = imposeRelpermScaling(model, 'SWL', .1, 'SWCR', .2, 'SWU', .9,...
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
state0     = initState(model.G, W , 350*barsa,[0.05,0.95]);

% Simulation tolerances
model.toleranceCNV = 1e-6;
%model.FacilityModel.toleranceWellRate = 1e-4;
%model.FacilityModel.toleranceWellMS  = 1e-4;

% Simulation schedule
% We use the first 70 time steps from the fine-scale simulation to train
% the network model. The remaining time steps will be used for prediction.
schedule = simpleSchedule(scheduleRef.step.val(1:70), 'W', W);

% Set up the problem structure
prob = struct('model', model, 'schedule', schedule, 'state0', state0);                               

% Plot the initial data
%figure
%plotCellData(model.G, state0.s(:,2))
%colorbar, view(0,0);axis equal tight;  daspect([1,0.1,0.1])


%% Specify the parameters to be matched
% We need to specify which of the parameters in the model we wish to match
% and any scaling and box limits that will applied to these in the BFGS
% optimization. These are specified through through the ModelParameter
% class from the optimization module.

% Parameter lumping
% In the GPSNet type of models, we only match a single pore volume and a
% single transmissibility for each network edge. We thus need a mapping
% between edges and cells/faces in the Cartesian grid to describe the
% corresponding lumping of parameters in this grid.
[cellEdgeNo, cellIx] = reorganizeIndices(ntwrkModel.Graph.Edges.Cell_Indices);
[faceEdgeNo, faceIx] = reorganizeIndices(ntwrkModel.Graph.Edges.Face_Indices);

nc =  model.G.cells.num;
nf =  numel(model.operators.T);

% Select configuration for sensitivity computations
config = {
    ...%name      include     scaling    boxlims    lumping      subset   relativeLimits
    'porevolume',       1,   'linear',       [],  cellEdgeNo,   cellIx,  [.001 2]
    'conntrans',        1,   'log',          [],          [],       [],  [.001 20]
    'transmissibility', 1,   'log'   ,       [],  faceEdgeNo,   faceIx,  [.001 20]
    'swl',              1,   'linear',       [],  ones(nc,1),       [],  [.5 2]
    'swcr',             1,   'linear',       [],  ones(nc,1),       [],  [.5 2]
    'swu',              1,   'linear',       [],  ones(nc,1),       [],  [.5 2]
    'sowcr',            1,   'linear',       [],  ones(nc,1),       [],  [.5 2]
    'krw',              1,   'linear',       [],  ones(nc,1),       [],  [.5 2]
    'kro',              1,   'linear',       [],  ones(nc,1),       [],  [.5 2]
    'sw',               1,   'linear',   [0 .3],  cellEdgeNo,   cellIx,  [0 10]};
    %'pressure'          1,   'linear',       [],  cellEdgeNo,   cellIx,  [.5 2]};
parameters = [];
for k = 1:size(config,1)
    if config{k, 2} == 0, continue, end
    parameters = addParameter(parameters, prob, ...
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
weighting =  {'WaterRateWeight',  day/4, ...
              'OilRateWeight',    day/4, ...
              'BHPWeight',        .5/barsa};
mismatchFn = @(model, states, schedule, states_ref, tt, tstep, state) ...
    matchObservedOW(model, states, schedule, states_ref,...
                    'computePartials', tt, 'tstep', tstep, weighting{:},...
                    'state',state,'from_states',false);

%% Evaluate the mismatch for the initial model
% To evaluate the initial mismatch, we first extract the network parameters
% from the initial specification, scale them to the unit interval, and pass
% them onto a library function from the optimization module that performs
% the forward simulation and evaluates the mismatch using mismatchFn.

% Extract and scale parameters
% Remember to overwrite pore volumes and transmissibilities if initial
% values for these have been computed by flow diagnostics
values = applyFunction(@(p)p.getParameterValue(prob), parameters);
values{2} =  5*values{2};
if any(strcmp(networkType,{'fd_preprocessor','fd_postprocessor'}))
    values{1} =  pv/cellsPerPath;
    values{3} =  TT;
end
for k = numel(values):-1:1    
    u{k} = parameters{k}.scale(values{k});
end
paramVec = vertcat(u{:});

% Perform a forward simulation and evaluate the mismatch
[misfitVal_0,~,wellSols_0,states_0] = ...
    evaluateMatch(paramVec,mismatchFn,prob,parameters, statesRef,'Gradient','none');

% Plot well curves for the reference and the initial model
plotWellSols({wellSolsRef, wellSols_0},...
    {scheduleRef.step.val, schedule.step.val}, ...
    'datasetnames', {'reference','initial'}, 'zoom', true, ...
    'field', 'qOs', 'SelectedWells', 9);


%% Model calibration
% Calibrate the model using the BFGS method
objh = @(p) evaluateMatch(p, mismatchFn, prob, parameters, statesRef);

[v, p_opt, history] = unitBoxBFGS(paramVec, objh,'objChangeTol',  1e-8, ...
    'maxIt',30, 'lbfgsStrategy', 'dynamic', 'lbfgsNum', 5);

%% Evaluate mismatch over the full simulation schedule
schedule = simpleSchedule(scheduleRef.step.val, 'W', W);
prob.schedule = schedule;

[~,~,wellSols_opt] = evaluateMatch(p_opt,   mismatchFn,prob,parameters,statesRef,'Gradient','none');
[~,~,wellSols_0]   = evaluateMatch(paramVec,mismatchFn,prob,parameters,statesRef,'Gradient','none');

%% Plot well curves
plotWellSols({wellSolsRef,wellSols_0,wellSols_opt}, ...
    {scheduleRef.step.val,schedule.step.val,schedule.step.val}, ...
    'datasetnames', {'reference','initial','optimized'}, 'zoom', true, ...
    'field', 'qOs', 'SelectedWells', 9);
legend('reference model','initial DD model','optimize DD model')                    