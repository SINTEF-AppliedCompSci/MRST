%% Calibration of a GPSNet model for SAIGUP using adjoints
% This example demonstrates how to set up a GPSNet type reduced network
% model based on injector-productor connections for a reservoir with
% realistic geology, but somewhat simplified flow physics and well
% placement. To this end, we use a realization from the SAIGUP project,
% which describes a shallow marine reservoir typical of the North Sea using
% a 40x120x20 corner-point grid with 78200 active cells.
%
% To calibrate the model, we use the BroydenFletcherGoldfarbShanno
% (BFGS) algorithm. This is an iterative line-search method that gradually
% improves an approximation to the Hessian matrix of the mismatch function,
% obtained only from adjoint gradients via a generalized secant method.
%
% This example was first introduced in MRST 2021b.
mrstModule add ad-core ad-blackoil deckformat diagnostics mrst-gui ...
    ad-props incomp optimization test-suite linearsolvers 

%% Setup 3D reference model
% Our reference model is the standard SAIGUP test case from the
% example-suite tutorial. This specifies an oil-water fluid model with a
% 5:1 viscosity ratio, quadratic relative permeabilities, and slight fluid
% compressibilities. Oil and water are initially separated by a sharp
% interface and the reservoir is produced using six vertical producers
% supported by eight vertical injectors around the perimeter of the model.
% If simulation results are not already stored for the model, we perform a
% full 3D simulation first.

predCase  = TestCase('saigup_wo');
predProbl = predCase.getPackedSimulationProblem();
%clearPackedSimulatorOutput(trueProb)
simulatePackedProblem(predProbl);

[predWellSols, predStates] = getPackedSimulatorOutput(predProbl);
predModel    = predCase.model;
predSchedule = predProbl.SimulatorSetup.schedule;
Wpred        = predSchedule.control.W;

% Random schedule
trainCase  = makeRandomTraining(predCase, 0.25, 0.05, false);
trainProbl = trainCase.getPackedSimulationProblem();
%clearPackedSimulatorOutput(trainProb)
simulatePackedProblem(trainProbl);

[trainWellSols, trainStates] = getPackedSimulatorOutput(trainProbl);
trainModel     = trainCase.model;
trainSchedule  = trainProbl.SimulatorSetup.schedule;
Wtrain         = trainSchedule.control.W;

%% Plot
predCase.plot(predStates,'step_index',numel(predStates));

plotWellSols({trainWellSols, predWellSols}, ...
    {trainSchedule.step.val, predSchedule.step.val},...
    'datasetnames',{'training','reference'}, ...
    'zoom', true, 'field', 'qWs', 'SelectedWells', 1:8);

%% Create the network
% We start by creating a network that connects injectors and producers.
% This network describes the possible 1D flow paths that each will have an
% associated pore volume and transmissibility.

% To create the network, we only need a single perforation from each well,
% which we somewhat arbitrarily pick to be the 7th from the top.
Wnw = Wpred;
for i = 1:numel(Wnw)
    Wnw(i).cells = Wnw(i).cells(7);
end

% Select type and create network 
% The module supplies different ways to do this: all to all, injectors to
% producers (and vice versa), connections specified manually by the user,
% connections determined by a flow diagnostics analysis of the 3D
% geological model, or connections based on flow diagnostics analysis of
% the fine-scale reference simulation.
networkType  = 'injectors_to_producers';
switch networkType
    case 'all_to_all'
        ntwrk =  Network(Wnw, predModel.G, 'type', networkType);
    case 'injectors_to_producers'
        ntwrk =  Network(Wnw, predModel.G, 'type', networkType,  ...
                         'injectors', 1:8, 'producers', 9:14);
    case 'fd_preprocessor'
        ntwrk =  Network(Wnw,predModel.G, 'type', networkType,  ...
                         'problem', predProbl, 'flow_filter',1*stb/day);
    case 'fd_postprocessor'
        ntwrk =  Network(Wnw,predModel.G, 'type', networkType,  ...
                         'problem', predProbl,   ...
                         'state_number',20,    ...
                         'flow_filter', 1*stb/day);
    otherwise
        error('\nType of network: %s is not implemented\n', networkType);             
end
                     
%% Create the GPSNet
% Each edge in the network is subgridded and mapped onto a row in a
% rectangular Cartesian gird having the same number of rows as the number
% of network edges. The fluid model is copied from the reference model.
gravity off
%S0 = mean(statesTrain{1}.s); S0=S0/sum(S0);
p0 = mean(trainStates{1}.pressure);
gpsNet = GPSNet(predModel, ntwrk, Wtrain, 'p0', p0);%, 'S0', S0);

% Reset the initial saturation, so that we have water in the left half of
% the model and oil in the right half. This will better reflect the
% situation with injectors completed in the water zone.
[I,~,~]=gridLogicalIndices(gpsNet.model.G);
gpsNet.state0.s(I<=5,:)=repmat([1 0], sum(I<=5),1);
figure, gpsNet.plotGrid(gpsNet.state0.s(:,1)); colorbar

%% Plot the GPSNet model: network and simulation grid
%fig1 = figure;
clf
subplot(2,2,1)
G = predCase.getVisualizationGrid();
plotCellData(G, log10(predCase.model.rock.perm(:,1)),'EdgeColor','none');
plotWell(G, predCase.schedule.control(1).W,'color','k','FontSize',10);
view(-85,65); 
axis tight off; set(gca,'Clipping',false)

subplot(2,2,3), ax=gca;
cmap=ntwrk.plotNetwork();
axis tight off; set(ax,'Clipping',false)
ax.Children(2).NodeFontSize=10;
view(-90,90);

subplot(2,2,4), ax=gca;
ntwrk.plotNetwork('circle'); 
axis equal tight off; set(ax,'Clipping',false)
ax.Children(1).NodeFontSize=10;

subplot(2,2,2), ax=gca;
gpsNet.plotGrid([]);
colormap(ax,min(cmap+.2,1));
for i=1:6, ax.Children(i).Position(1) = 3100; end
for i=7:14, ax.Children(i).Position(1)= -150; end
set(ax.Children(1:14),'FontSize',10);

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

% Select configuration for sensitivity computations
config = {
    ...%name      include     scaling    boxlims    lumping      subset   relativeLimits
    'porevolume',       1,   'linear',       [],  cellEdgeNo,   cellIx,  [.001 2]
    'conntrans',        1,   'log',          [],          [],       [],  [.001 20]
    'transmissibility', 1,   'log'   ,       [],  faceEdgeNo,   faceIx,  [.001 20]
    'sw',               1,   'linear',    [0,1],          [],       [],  []};
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
weighting =  {'WaterRateWeight',  day/2000, ...
              'OilRateWeight',    day/1000, ...
              'BHPWeight',        1/(50*barsa)};
mismatchFn = @(model, states, schedule, states_ref, tt, tstep, state) ...
    matchObservedOW(model, states, schedule, states_ref,...
                    'computePartials', tt, 'tstep', tstep, weighting{:},...
                    'state',state,'from_states',false);

%% Evaluate the initial misfit in training data and prediction data
% To define the initial model, we extract the network parameters from the
% initial specification, scale them to the unit interval, and organize them
% into a vector that can be passed onto mismatchFn
pinit = gpsNet.getScaledParameterVector(trainSetup, trainPrms, 'connscale', 5);
 
% Evaluate the initial mismatch for the training data
[misfitT0,~,wellSolT0] = ...
    evaluateMatch(pinit,mismatchFn,trainSetup,trainPrms,trainStates,'Gradient','none');
[misfitP0,~,wellSolP0] = ...
    evaluateMatch(pinit,mismatchFn,predSetup,predPrms,predStates,'Gradient','none');

%% Model calibration
% Calibrate the model using the BFGS method. This is a computationally
% expensive operation that may run for several hours if you choose a large
% number of iterations. Here, we therefore only apply 10 iterations. 
objh = @(p)evaluateMatch(p,mismatchFn,trainSetup,trainPrms,trainStates);
[v, popt, history] = unitBoxBFGS(pinit, objh, 'objChangeTol', 1e-8, ...
    'gradTol', 1e-5, 'maxIt', 10, 'lbfgsStrategy', 'dynamic', ...
    'lbfgsNum', 5, 'outputHessian', true, 'logPlot', true);
    
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
    'field', 'qWs', 'SelectedWells', [9,14]);
set(fh, 'name','SAIGUP: GPSNet training')


%% Plot well responses for prediction case
predSteps = predSchedule.step.val;
fh = plotWellSols({predWellSols,wellSolP0, wellSolP}, ...
    {predSteps, predSteps, predSteps}, ...
    'datasetnames', {'reference','initial','predicted'}, 'zoom', true, ...
    'field', 'qWs', 'SelectedWells', [9,14]);
set(fh, 'name','SAIGUP: GPSNet prediction')


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
