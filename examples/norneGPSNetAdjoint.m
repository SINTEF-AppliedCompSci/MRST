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
% This example was first introduced in MRST 2021b.
mrstModule add ad-core ad-blackoil deckformat diagnostics...
               mrst-gui ad-props incomp optimization...
               network-models test-suite linearsolvers 

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
predCase  = TestCase('norne_simple_wo');
predProbl = predCase.getPackedSimulationProblem();
%clearPackedSimulatorOutput(trueProb)
simulatePackedProblem(predProbl);

[predWellSols, predStates] = getPackedSimulatorOutput(predProbl);
predModel     = predCase.model;
predSchedule  = predProbl.SimulatorSetup.schedule;
Wpred         = predSchedule.control.W;

% Random schedule
trainCase  = makeRandomTraining(predCase, 0.25, 0.05, false);
trainProbl = trainCase.getPackedSimulationProblem();
%clearPackedSimulatorOutput(trainProb)
simulatePackedProblem(trainProbl);

[trainWellSols, trainStates] = getPackedSimulatorOutput(trainProbl);
trainModel     = trainCase.model;
trainSchedule  = trainProbl.SimulatorSetup.schedule;
Wtrain         = trainSchedule.control.W;

% Plot
predCase.plot(predStates,'step_index',numel(predStates))

plotWellSols({trainWellSols, predWellSols}, ...
    {trainSchedule.step.val, predSchedule.step.val},...
    'datasetnames',{'training','reference'}, ...
    'zoom', true, 'field', 'qWs', 'SelectedWells', 1:6);


%% Create the network
% We start by creating a network that connects injectors and producers.
% This network describes the possible 1D flow paths that each will have an
% associated pore volume and transmissibility.

% To create the network, we only need a single perforation from each well,
% which we somewhat arbitrarily pick to be the middle perforation.
Wnw = Wpred;
for i = 1:numel(Wnw)
    Wnw(i).cells = Wnw(i).cells(round(numel(Wnw(i).cells)/2));
end

% Select and create network
% The module offers different ways to specify the network topology. To
% choose different ones, you change the networkType variable
networkType = 'injectors_to_producers';
switch networkType
    case 'all_to_all'
        ntwrk =  Network(Wnw, predModel.G, 'type', networkType);
    case 'injectors_to_producers'
        ntwrk =  Network(Wnw, predModel.G, 'type', networkType, ...
                         'injectors', 1:6, 'producers', 7:11);
    case 'user_defined_edges'
        edges = [rldecode((1:6)',5*ones(6,1)), repmat((7:11)',6,1)];
        edges = sortrows([edges; edges],'ascend');
        ntwrk =  Network(Wnw, predModel.G, 'type', networkType, ...
                         'edges',edges);
    case 'fd_preprocessor'
        ntwrk =  Network(Wnw, predModel.G, 'type', networkType, ...
                         'problem', predProbl, 'flow_filter', 1*stb/day);
    case 'fd_postprocessor'
        ntwrk =  Network(Wnw, predModel.G, 'type', networkType, ...
                         'problem', predProbl,  ...
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
gpsNet = GPSNet(predModel, ntwrk, Wtrain);

%% Plot the GPSNet model: network and simulation grid
fig1 = figure;
subplot(2,2,1)
G = predCase.getVisualizationGrid();
plotCellData(G, predCase.model.rock.poro,'EdgeColor','none');
plotWell(G, predCase.schedule.control(1).W,'color','k','FontSize',10);
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
trainSetup = gpsNetSimulationSetup(gpsNet, trainSchedule);
predSetup  = gpsNetSimulationSetup(gpsNet, predSchedule);

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
trainPrms = []; predPrms = [];
for k = 1:size(config,1)
    if config{k, 2} == 0, continue, end
    trainPrms = addParameter(trainPrms, trainSetup, ...
        'name',    config{k,1}, 'scaling', config{k,3}, ...
        'boxLims', config{k,4}, 'lumping', config{k,5}, ...
        'subset',  config{k,6}, 'relativeLimits',config{k,7});
    predPrms = addParameter(predPrms, predSetup, ...
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
pinit = gpsNet.getScaledParameterVector(trainSetup, trainPrms, 'connscale', 0.5);
 
% Evaluate the initial mismatch for the training data
[misfitT0,~,wellSolT0] = ...
    evaluateMatch(pinit,mismatchFn,trainSetup,trainPrms,trainStates,'Gradient','none');
[misfitP0,~,wellSolP0] = ...
    evaluateMatch(pinit,mismatchFn,predSetup,predPrms,predStates,'Gradient','none');

%% Model calibration
% Calibrate the model using the BFGS method. This is a computationally
% expensive operation that may run for several hours if you choose a large
% number of iterations. Here, we therefore only apply 10 iterations. To get
% a good match, you should increase this number to 100+.  
objh = @(p)evaluateMatch(p,mismatchFn,trainSetup,trainPrms,trainStates);
[v, popt, history] = unitBoxBFGS(pinit, objh, 'objChangeTol', 1e-8, 'gradTol', 1e-5, ...
    'maxIt', 10, 'lbfgsStrategy', 'dynamic', 'lbfgsNum', 5, ...
    'outputHessian', true, 'logPlot', true);

% If you are not happy with the match you have obtained so far, you can
% continue iterating if you supply the history as an optional parameter as
% follows:
%[v, popt, history] = unitBoxBFGS(popt, objh, 'objChangeTol', 1e-8, 'gradTol', 1e-5, ...
%    'maxIt', 10, 'lbfgsStrategy', 'dynamic', 'lbfgsNum', 5, ...
%    'outputHessian', true, 'logPlot', true, 'history', history);

%% Evaluate mismatch over the full simulation schedule 
[misfitT,~,wellSolT] = ...
    evaluateMatch(popt,mismatchFn,trainSetup,trainPrms,trainStates,'Gradient','none');
[misfitP, ~,wellSolP] = ...
    evaluateMatch(popt,mismatchFn,predSetup,predPrms,predStates,'Gradient','none');

%% Plot well responses for training data
trainSteps = trainSchedule.step.val;
fh = plotWellSols({trainWellSols,wellSolT0,wellSolT}, ...
    {trainSteps, trainSteps, trainSteps}, ...
    'datasetnames', {'train','init','match'}, 'zoom', true, ...
    'field', 'qWs', 'SelectedWells', 7:8);
set(fh, 'name','Norne: GPSNet training')


%% Plot well responses for prediction case
predSteps = predSchedule.step.val;
fh = plotWellSols({predWellSols,wellSolP0, wellSolP}, ...
    {predSteps, predSteps, predSteps}, ...
    'datasetnames', {'reference','initial','predicted'}, 'zoom', true, ...
    'field', 'qOs', 'SelectedWells', 7:8);
set(fh, 'name','Norne: GPSNet prediction')

%% Copyright Notice
%
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
