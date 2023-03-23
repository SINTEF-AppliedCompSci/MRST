%% Norne: calibration of a GPSNet model
% This example demonstrates how to set up a GPSNet type reduced network
% model for a geostatistical realization of the Norne field model with
% simplified flow physics and well placement.
%
% To calibrate the model, we specify a training simulation with oscillatory
% well controls and use the BroydenFletcherGoldfarbShanno (BFGS)
% algorithm. This is an iterative line-search method that gradually
% improves an approximation to the Hessian matrix of the mismatch function,
% obtained only from adjoint gradients via a generalized secant method.
%
% This example was first introduced in MRST 2022a.
mrstModule add ad-core ad-blackoil deckformat diagnostics...
               mrst-gui ad-props incomp optimization...
               network-models test-suite linearsolvers 
         
%% Setup 3D reference model
% The reference model is a single stochastic realization of the OLYMPUS
% field model (the first ensemble member). Because MRST does not support
% the TUNING keyword in ECLIPSE input decks, we set the schedule manually.

% To train the GPSNet model, we set up a schedule that oscillates randomly
% around the prescribed pressure controls of the reference solution. This
% is done to generate a more general model that should provide accurate
% predictions in a sufficiently small region around the specified controls.

% True setup, which we seek to reproduce
predCase  = MRSTExample('olympus_field_wo');
predProbl = predCase.getPackedSimulationProblem();
%clearPackedSimulatorOutput(problem)
simulatePackedProblem(predProbl);

[wellSolPred, statesPred] = getPackedSimulatorOutput(predProbl);
predModel    = predCase.model;
predSchedule = predProbl.SimulatorSetup.schedule;
Wpred        = predSchedule.control.W;

% Random schedule
trainCase  = makeRandomTraining(predCase, 0.0, 0.05, false);
trainProbl = trainCase.getPackedSimulationProblem();
%clearPackedSimulatorOutput(trainProb)
simulatePackedProblem(trainProbl);

[wellSolTrain, statesTrain] = getPackedSimulatorOutput(trainProbl);
trainModel     = trainCase.model;
trainSchedule  = trainProbl.SimulatorSetup.schedule;
Wtrain         = trainSchedule.control.W;

% Plot
%trueEx.plot(statesTrue,'step_index',numel(statesTrue))

plotWellSols({wellSolTrain, wellSolPred}, ...
    {trainSchedule.step.val, predSchedule.step.val},...
    'datasetnames',{'training','reference'}, ...
    'zoom', true, 'field', 'qWs', 'SelectedWells', 1:7);


%% Create the network
% We start by creating a network that connects injectors and producers.
% This network describes the possible 1D flow paths that each will have an
% associated pore volume and transmissibility. Because the original model
% has multiple fluid regions, we cannot use flow diagnostics preprocessing,
% since the incompressible flow solvers this relies on do not yet handle
% such fluid object. 
Wnw = Wpred;
for i = 1:numel(Wnw)
    Wnw(i).cells = Wnw(i).cells(round(numel(Wnw(i).cells)/2));
end
networkType = 'injectors_to_producers';
switch networkType
    case 'all_to_all'
        ntwrk =  Network(Wnw, predModel.G, 'type', networkType);
    case 'injectors_to_producers'
        ntwrk =  Network(Wnw, predModel.G, 'type', networkType, ...
                         'injectors', 1:7, 'producers', 8:18);
    case 'fd_postprocessor'
        ntwrk =  Network(Wnw, predModel.G, 'type', networkType, ...
                         'problem', predProbl,                   ...
                         'flow_filter', 1*stb/day,             ...
                         'state_number', numel(statesPred));
    otherwise
        error('\nNetwork of type %s is not implemented\n', networkType);           
end

%% Create the GPSNet
% Each edge in the network is subgridded and mapped onto a row in a
% rectangular Cartesian grid having the same number of rows as the number
% of network edges. By default, the fluid model is copied from the
% fine-scale model. The GPSNet is not yet programmed to setup multiple
% fluid regions and will hence only pick the fluid model from the first
% region. This may be ok.
gravity off
gpsNet = GPSNet(predModel, ntwrk, Wtrain);

%{
% If this is not ok, we can replace the fluid model with a custom-made
% fluid model having quadratic relative permeabilities, which can then be
% further calibrated locally for each connection by adjusting the endpoint
% scaling.
fluid = initSimpleADIFluid('phases', 'WO',... 
                           'mu' , [.3, 3]*centi*poise,...
                           'rho', [1014, 859]*kilogram/meter^3, ...
                           'n', [2 2]);
fluid.krPts  = struct('w', [0 0 1 1], 'ow', [0 0 1 1]);
gpsNet = NetworkModel(modelTrue, ntwrk, WTrain, 'fluid', fluid, ...
    'scaling', {'SWL', .1, 'SWCR', .2, 'SWU', .9, ...
                'SOWCR', .1, 'KRW', .9, 'KRO', .8});
%}

%% Plot the GPSNet model: network and simulation grid
fig1 = figure;
subplot(2,2,1)
G = predCase.getVisualizationGrid();
plotCellData(G, predCase.model.rock.poro,'EdgeColor','none');
plotWell(G, predCase.schedule.control(1).W,'color','k','FontSize',10);
view(20,60); 
axis tight off; set(gca,'Clipping',false); zoom(1.2);

subplot(2,2,3), ax=gca;
cmap=ntwrk.plotNetwork('onGrid',false);
axis tight off; set(ax,'Clipping',false); zoom(1.3); 

subplot(2,2,4), ax=gca;
ntwrk.plotNetwork('circle'); 
axis equal tight off; set(ax,'Clipping',false); zoom(1.3); 

subplot(2,2,2), ax=gca;
gpsNet.plotGrid([]);
colormap(ax,min(cmap+.2,1));
for i=1:11, ax.Children(i).Position(1:2) =  [2830,-10]; end
for i=16:18, ax.Children(i).Position(1:2) = [20,-10]; end
set(ax.Children(1:11),'FontSize',8);
ax.XLim = [-125 3300]; ax.ZLim = [0 200];
ax.Position([2 4]) = ax.Position([2 4])+[-.05 .1];

%% Specify training schedule and parameters to be matched
% We need to specify which of the parameters in the model we wish to match
% and any scaling and box limits that will applied to these in the BFGS
% optimization. These are specified through through the ModelParameter
% class from the optimization module. For the training we use the random
% schedule
trainSetup = gpsNetSimulationSetup(gpsNet, trainSchedule);
predSetup  = gpsNetSimulationSetup(gpsNet, predSchedule);

% Parameter lumping
% In the GPSNet type of models, we only match a single pore volume and a
% single transmissibility for each network edge. We thus need a mapping
% between edges and all associated cells/faces in the Cartesian grid to
% describe the corresponding lumping of parameters in this grid.
[cellEdgeNo, cellIx] = gpsNet.getMapping('cells');
[faceEdgeNo, faceIx] = gpsNet.getMapping('faces');

nc =  gpsNet.model.G.cells.num;
nf =  numel(gpsNet.model.operators.T);

% select configuration for sensitivity computations
config = {
    ...%name      include    scaling     lumping   subset     boxlims  relativeLimits
    'porevolume',       1,  'linear',  cellEdgeNo,  cellIx,        [],  [1e-2 10]
    'conntrans',        1,  'log',             [],      [],        [],  [1e-2 1e4]
    'transmissibility', 1,  'log'   ,  faceEdgeNo,  faceIx,        [],  [1e-2 1e3]
    'swl',              1,  'linear',  cellEdgeNo,  cellIx,    [0,.4],  []
    'swcr',             1,  'linear',  cellEdgeNo,  cellIx,    [0,.4],  []
    'swu',              1,  'linear',  cellEdgeNo,  cellIx,    [.8,1],  []
    'sowcr',            1,  'linear',  cellEdgeNo,  cellIx,    [0,.4],  []
    'krw',              1,  'linear',  cellEdgeNo,  cellIx,  [.5,1.5],  []
    'kro',              1,  'linear',  cellEdgeNo,  cellIx,  [.5,1.5],  []
    'sw',               1,  'linear',          [],      [],     [0,1],  []};
trainPrms = []; predPrms = [];
for k = 1:size(config,1)
    if config{k, 2} == 0, continue, end
    trainPrms = addParameter(trainPrms, trainSetup, ...
        'name',    config{k,1}, 'scaling',        config{k,3}, ...
        'lumping', config{k,4}, 'subset',         config{k,5}, ...
        'boxLims', config{k,6}, 'relativeLimits', config{k,7});
    predPrms  = addParameter(predPrms,  predSetup, ...
        'name',    config{k,1}, 'scaling',        config{k,3}, ...
        'lumping', config{k,4}, 'subset',         config{k,5}, ...
        'boxLims', config{k,6}, 'relativeLimits', config{k,7});
end

%% Define the mismatch function
% The mismatch function is defined as a function handle to a library
% function from the optimization module that computes the mismatch between
% a given simulation and a reference state. For an oil-water system, the
% match is computed based on three quantities (water/oil rate and bhp) and
% these must be given an associated weight. Because the wells are
% controlled by bhp, we only match rates here.
weighting =  {'WaterRateWeight',  day/2000, ...
              'OilRateWeight',    day/900};    
mismatchFn = @(model, states, schedule, states_ref, tt, tstep, state) ...
    matchObservedOW(model, states, schedule, states_ref,...
                    'computePartials', tt, 'tstep', tstep, weighting{:},...
                    'state', state, 'from_states', false);
          
%% Evaluate the initial misfit in training data and prediction data
% To define the initial model, we extract the network parameters from the
% initial specification, scale them to the unit interval, and organize them
% into a vector that can be passed onto mismatchFn. Because the fine model
% has a nonuniform initial saturation distribution, it is important to
% initialize the initial saturation in the GPSNet randomly. 
pinit = gpsNet.getScaledParameterVector(trainSetup, trainPrms, 'randSw', true);
 
% Evaluate the initial mismatch for the training data
[misfitT0,~,wellSolT0] = ...
    evaluateMatch(pinit,mismatchFn,trainSetup,trainPrms,statesTrain,'Gradient','none');
[misfitP0,~,wellSolP0] = ...
    evaluateMatch(pinit,mismatchFn,predSetup, predPrms, statesPred, 'Gradient','none');
 
%% Model calibration with BFGS
objh = @(p)evaluateMatch(p,mismatchFn,trainSetup,trainPrms,statesTrain);
[v, popt, history] = unitBoxBFGS(pinit, objh,'objChangeTol', 1e-8, ...
    'gradTol', 1e-5, 'maxIt', 10, 'lbfgsStrategy', 'dynamic', ...
    'lbfgsNum', 5, 'outputHessian', true, 'logPlot', true);

%% Evaluate mismatch over the full simulation schedule 
[misfitT,~,wellSolT] = ...
    evaluateMatch(popt,mismatchFn,trainSetup,trainPrms,statesTrain,'Gradient','none');
[misfitP, ~,wellSolP] = ...
    evaluateMatch(popt,mismatchFn,predSetup,predPrms,statesPred,'Gradient','none');

%% Plot well responses for training data
trainSteps = trainSchedule.step.val;
fh = plotWellSols({wellSolsTrain,wellSolT0,wellSolT}, ...
    {trainSteps, trainSteps, trainSteps}, ...
    'datasetnames', {'train','init','match'}, 'zoom', true, ...
    'field', 'qWs', 'SelectedWells', 1);
set(fh, 'name','OLYMPUS: GPSNet training')


%% Plot well responses for prediction case
predSteps = predSchedule.step.val;
fh = plotWellSols({wellSolPred,wellSolP0, wellSolP}, ...
    {predSteps, predSteps, predSteps}, ...
    'datasetnames', {'reference','initial','predicted'}, 'zoom', true, ...
    'field', 'qWs', 'SelectedWells', 1);
set(fh, 'name','OLYMPUS: GPSNet prediction')

%% Copyright Notice
%
% <html>
% <p><font size="-1">
% Copyright 2009-2023 SINTEF Digital, Mathematics & Cybernetics.
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
