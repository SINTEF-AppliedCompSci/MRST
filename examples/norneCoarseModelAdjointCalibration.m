%% Norne: calibration of a network model formed from a coarse grid
% This example demonstrates how to set up a data driven model with 3D type
% graph topology and corresponding parameters initialized based on a coarse
% partition of the full corner-point grid of the Norne field model.
%
% To calibrate the model, we use the Broyden–Fletcher–Goldfarb–Shanno
% (BFGS) algorithm. This is an iterative line-search method that gradually
% improves an approximation to the Hessian matrix of the mismatch function,
% obtained only from adjoint gradients via a generalized secant method.
mrstModule add ad-core ad-blackoil deckformat ...
               agglom upscaling coarsegrid...
               mrst-gui ad-props incomp optimization...
               network-models example-suite linearsolvers 

%% Setup 3D reference model
% The reference model is a single stochastic realization of the Norne field
% model taken from the example-suit module. Compared with the real field,
% this simulation case has simpler fluid description (an oil-water model)
% and an idealized field development plan consisting of a simple pattern of
% eleven vertical wells that run under constant bhp or rate controls.

example = MRSTExample('norne_simple_wo');

% Modify the setup to use shorter time steps
problem  = example.getPackedSimulationProblem();
schedule = example.schedule;
ix  = repmat(1:numel(schedule.step.val), [4 1]);
schedule.step.val     = schedule.step.val(ix(:))/4;
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

%% Coarse-scale model
% We make a coarse grid defined by a uniform 6 x 8 x 1 partition and
% perform a simple upscaling to obtain a coarse model
cModel = upscaleModelTPFA(modelRef, partitionUI(modelRef.G,[6 8 1]));

% We want to include rel-perm scaling as tunabale parameters, so include
% these for the coarse model. These parameters have no effect for the
% initial coarse model (they are set equal to the ones given by the
% rel-perm curves).
pts = cModel.fluid.krPts;
scaling = {'SWL',   pts.w(1), 'SWCR', pts.w(2), 'SWU', pts.w(3), ...
           'SOWCR', pts.o(2), 'KRW',  pts.w(4), 'KRO', pts.o(4)};
cModel = imposeRelpermScaling(cModel, scaling{:});
cModel.toleranceCNV = 1e-6;  % tighter tolerance to improve gradient accuracy

%% Plot comparison of the coarse and fine model
figure('position',[100 100 1000 400])
axes('position',[.02 .05 .48 .9]);
plotCellData(modelRef.G,modelRef.rock.poro,'EdgeAlpha',.1); 
view(174,60);
title('Fine-scale model porosity')
plotWell(modelRef.G,Wref,'Color','k'); axis off tight
mrstColorbar(modelRef.rock.poro,'South'); cax = caxis;

axes('position',[.5 .05 .48 .9]);
plotCellData(cModel.G, cModel.rock.poro,'EdgeColor','none');
title('Coarse-scale model porosity')
plotFaces(cModel.G, boundaryFaces(cModel.G), 'EdgeColor', [0.4 0.4 0.4], ...
         'EdgeAlpha',.5, 'FaceColor', 'none');
view(174,60); caxis(cax);
plotWell(modelRef.G, Wref, 'Color', 'k'); axis off tight
mrstColorbar(cModel.rock.poro,'South');

%% Simulate initial upscaled coarse model for full time
cState0   = upscaleState(cModel, modelRef, example.state0);
cSchedule = upscaleSchedule(cModel, scheduleRef);
[cWellSols, cStates] = simulateScheduleAD(cState0, cModel, cSchedule);
plotWellSols({wellSolsRef, cWellSols}, ...
    {scheduleRef.step.val, cSchedule.step.val},...
    'datasetnames',{'fine scale model','initial upscaled model'}, ...
    'zoom', true, 'field', 'qOs', 'SelectedWells', 7);
drawnow

%% Specify training schedule and parameters to be matched
% We use the first half of the given data for training. In this setup, we
% use all pore volumes, transmissibilities, and well connections in the
% coarse grid as calibration parameters.
trainSteps = 1:round(numel(scheduleRef.step.val)/2);
timeSteps  = scheduleRef.step.val(trainSteps); 
trainSched = upscaleSchedule(cModel, simpleSchedule(timeSteps, 'W', Wref));
trainProbl = struct('model', cModel, 'schedule', trainSched, 'state0', cState0);

cellIx =  ones(cModel.G.cells.num,1);
config = {
    ...%name      include     scaling    boxlims  lumping   subset  relativeLimits
    'porevolume',       1,   'linear',       [],    [],      [],      [.001 4]
    'conntrans',        1,   'log',          [],    [],      [],      [.001 100]
    'transmissibility', 1,   'log'   ,       [],    [],      [],      [.001 100]
    'swl',              1,   'linear',       [],  cellIx,    [],       [.5 2]
    'swcr',             1,   'linear',       [],  cellIx,    [],       [.5 2]
    'swu',              1,   'linear',       [],  cellIx,    [],       [.5 2]
    'sowcr',            1,   'linear',       [],  cellIx,    [],       [.5 2]
    'krw',              1,   'linear',       [],  cellIx,    [],       [.5 2]
    'kro',              1,   'linear',       [],  cellIx,    [],       [.5 2]
    'sw',               1,   'linear',   [0 .6],    [],      [],       [0 10]
    'pressure'          1,   'linear',       [],    [],      [],       [.1 4]};
parameters = [];
for k = 1:size(config,1)
    if config{k, 2} == 0, continue, end
    parameters = addParameter(parameters, trainProbl, ...
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
weighting  = {'WaterRateWeight',  day/500, ...
              'OilRateWeight',    day/500, ...
              'BHPWeight',        1/(50*barsa)};
mismatchFn = @(model, states, schedule, states_ref, tt, tstep, state) ...
    matchObservedOW(model, states, schedule, states_ref,...
                   'computePartials', tt, 'tstep', tstep, weighting{:},...
                   'state', state, 'from_states', false);
       
%% Evaluate the mismatch for the initial model
% We extract and scale the initial parameter set and use this to evaluate
% the forward model and then measure the mismatch
%[misfitVal0,~,wellSols0,states0] = ...
%    evaluateMatch(pvec,obj,trainProbl,parameters, statesRef,'Gradient','none');          
 
% Plot solution
%plotWellSols({wellSolsRef,wellSols0},{scheduleRef.step.val,trainProbl.schedule.step.val})


%% Model calibration
pvec = getScaledParameterVector(trainProbl, parameters);
objh = @(p) evaluateMatch(p,mismatchFn,trainProbl,parameters,statesRef);

[v, p_opt, history] = unitBoxBFGS(pvec, objh, 'objChangeTol', 1e-8, ...
    'maxIt',30, 'lbfgsStrategy', 'dynamic', 'lbfgsNum', 5);

%% Evaluate mismatch over the full simulation schedule
prob = trainProbl;
prob.schedule = simpleSchedule(scheduleRef.step.val, 'W', prob.schedule.control.W);
 
[~,~,wellSols_opt] = evaluateMatch(p_opt,mismatchFn,prob,parameters,statesRef,'Gradient','none');
[~,~,wellSols0]    = evaluateMatch(pvec, mismatchFn,prob,parameters,statesRef,'Gradient','none');

%% Plot well curves
fh = plotWellSols({wellSolsRef,wellSols0,wellSols_opt}, ...
    {scheduleRef.step.val,prob.schedule.step.val,prob.schedule.step.val}, ...
    'datasetnames', {'reference','initial','optimized'}, 'zoom', true, ...
    'field', 'qOs', 'SelectedWells', 7);
set(fh, 'name','Norne')
legend('reference model','initial DD model','optimize DD model')

%% Plot the matched porosities
% In a network model, the calibrated parameters cannot generally be
% interpreted as in their original physical sense, e.g., because we do not
% make any attempts to preserve any geological realism. We nonetheless plot
% the calibrated porosities for comparison. 
figure
plotCellData(cModel.G, p_opt(1:cModel.G.cells.num),'EdgeColor','none');
plotFaces(cModel.G, boundaryFaces(cModel.G), 'EdgeColor', [0.4 0.4 0.4], ...
         'EdgeAlpha',.5, 'FaceColor', 'none');
view(174,60); caxis(cax);
plotWell(modelRef.G, Wref, 'Color', 'k'); axis off tight
mrstColorbar(p_opt(1:cModel.G.cells.num),'South');
                                       