%% Calibration of a GPSNet model for Egg using adjoints
% This example demonstrates how to set up a GPSNet type reduced network for
% the Egg model, a conceptual model of a channelized reservoir. For details
% about the model, see Jansen, J. D., et al. "The egg model -- a geological
% ensemble for reservoir simulation." Geoscience Data Journal 1.2 (2014):
% 192-195. The network topology can be determined in different ways. The
% default setup of the example is to use flow diagnostic postprocessing of
% the reference simulation to determine nonzero interwell connections.
%
% For the calibration, we use the limited-memory
% BroydenFletcherGoldfarbShanno (L-BFGS) algorithm. This is an iterative
% line-search method that gradually improves an approximation to the
% Hessian matrix of the mismatch function, obtained only from adjoint
% gradients via a generalized secant method.
%
% This example was first introduced in MRST 2021b.
mrstModule add ad-core ad-blackoil deckformat diagnostics mrst-gui ...
    ad-props incomp optimization network-models test-suite linearsolvers

%% Setup 3D reference model
% The Egg data set consists of a base case and an ensemble of one hundred
% different permeability realizations. (If you do not have the data
% available, you have to download it manually from the internet; see
% mrstDatasetGUI for details.) Herein, we will use realization number 0,
% the base case, but you can also run the example with any number between 1
% and 100.

% True schedule, which we seek to reproduce
predCase  = TestCase('egg_wo');
predProbl = predCase.getPackedSimulationProblem();
%clearPackedSimulatorOutput(trueCase)
simulatePackedProblem(predProbl);

[predWellSols, predStates] = getPackedSimulatorOutput(predProbl);
predModel     = predCase.model;
predSchedule  = predProbl.SimulatorSetup.schedule;
Wpred         = predSchedule.control.W;

% Plot
predCase.plot(predStates,'step_index',numel(predStates));


%% Random schedule
% The basic schedule produces the reservoir with bhp controls that are very
% close to the reservoir pressure. To avoid prescribing conditions that
% result in injection from the producers, we introduce a bhp perturbation
% that is nonsymmetric around the base case
trainCase  = makeRandomTraining(predCase, 0.25, @(x,y) y - 5*(x-.2)*barsa, false);
trainProbl = trainCase.getPackedSimulationProblem();
%clearPackedSimulatorOutput(trainCase)
simulatePackedProblem(trainProbl);

[trainWellSols, trainStates] = getPackedSimulatorOutput(trainProbl);
trainModel     = trainCase.model;
predSchedule  = trainProbl.SimulatorSetup.schedule;
Wtrain         = predSchedule.control.W;

plotWellSols({trainWellSols, predWellSols}, ...
    {predSchedule.step.val, predSchedule.step.val},...
    'datasetnames',{'training','reference'}, ...
    'zoom', true, 'field', 'qWs', 'SelectedWells', [1 9]);

%% Create the network
% We start by creating a network that connects injectors and producers.
% This network describes the possible 1D flow paths that each will have an
% associated pore volume and transmissibility.

% To create the network, we only need a single perforation from each well,
% which we somewhat arbitrarily pick from the top.
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
networkType = 'injectors_to_producers';
switch networkType
    case 'all_to_all'
        ntwrk =  Network(Wnw, predModel.G, 'type', networkType);
    case 'injectors_to_producers'
        ntwrk =  Network(Wnw, predModel.G, 'type', networkType, ...
                         'injectors', 1:8, 'producers', 9:12);
    case 'user_defined_edges'
        edges = [1 9;2 9;2 10;3 9;3 11;
                 4 9; 4 10; 4 11; 4 12;
                 5 10; 5 12; 6 11; 7 11; 7 12; 8 12];
        ntwrk =  Network(Wnw, predModel.G, 'type', networkType, ...
                         'edges',edges);
    case 'fd_preprocessor'
         ntwrk = Network(Wnw, predModel.G, 'type', networkType, ...
                         'problem', predProbl,                   ...
                         'flow_filter',1*stb/day);
    case 'fd_postprocessor'
        ntwrk = Network(Wnw, predModel.G, 'type', networkType, ...
                         'problem', predProbl,                 ...
                         'state_number',40,                   ...
                         'flow_filter', 1*stb/day);
    otherwise
        error('\nNetwork of type %s is not implemented\n', networkType);
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
plotCellData(G, log10(predCase.model.rock.perm(:,1)),'EdgeColor','none');
plotWell(G, predCase.schedule.control(1).W,'color','k','FontSize',8);
view(85,65); 
axis tight off; set(gca,'Clipping',false)

subplot(2,2,3), ax=gca;
cmap=ntwrk.plotNetwork();
axis tight off; set(ax,'Clipping',false)
ax.Children(2).NodeFontSize=8;

subplot(2,2,4), ax=gca;
ntwrk.plotNetwork('circle'); 
axis equal tight off; set(ax,'Clipping',false)
ax.Children(1).NodeFontSize=8;

subplot(2,2,2), ax=gca;
gpsNet.plotGrid([]);
colormap(ax,min(cmap+.2,1));
for i=1:4, ax.Children(i).Position(1:2) =  [495,50]; end
for i=5:12, ax.Children(i).Position(1:2) = [-60,50]; end
set(ax.Children(1:12),'FontSize',8);

%% Specify training schedule and parameters to be matched
% We need to specify which of the parameters in the model we wish to match
% and any scaling and box limits that will applied to these in the BFGS
% optimization. These are specified through through the ModelParameter
% class from the optimization module. For the training we use the random
% schedule
trainSetup = gpsNetSimulationSetup(gpsNet, predSchedule);
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
    'porevolume',       1,   'linear',       [],  cellEdgeNo,   cellIx,   [.01  5]
    'conntrans',        1,   'log',          [],          [],       [],   [.001 20]
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
weighting =  {'WaterRateWeight',  day/150, ...
              'OilRateWeight',    day/150, ...
              'BHPWeight',        1/(40*barsa)};   
mismatchFn = @(model, states, schedule, states_ref, tt, tstep, state) ...
    matchObservedOW(model, states, schedule, states_ref,...
                    'computePartials', tt, 'tstep', tstep, weighting{:},...
                    'state',state,'from_states',false);

%% Evaluate the initial misfit in training data and prediction data
% To define the initial model, we extract the network parameters from the
% initial specification, scale them to the unit interval, and organize them
% into a vector that can be passed onto mismatchFn
pinit = gpsNet.getScaledParameterVector(trainSetup, trainPrms, 'connscale', 7);
 
% Evaluate the initial mismatch for the training data
[misfitT0,~,wellSolT0] = ...
    evaluateMatch(pinit,mismatchFn,trainSetup,trainPrms,trainStates,'Gradient','none');
[misfitP0,~,wellSolP0] = ...
    evaluateMatch(pinit,mismatchFn,predSetup,predPrms,predStates,'Gradient','none');


%% Model calibration
% Calibrate the model using the BFGS method. This is a computationally
% expensive operation that may run for several hours if you choose a large
% number of iterations. Here, we therefore only apply 25 iterations, which
% may take up to ten minutes to run, depending on your computer.
objh = @(p)evaluateMatch(p,mismatchFn,trainSetup,trainPrms,trainStates);
[v, popt, history] = unitBoxBFGS(pinit, objh, 'objChangeTol', 1e-8, ...
    'gradTol', 1e-5, 'maxIt', 25, 'lbfgsStrategy', 'dynamic', ...
    'lbfgsNum', 5, 'outputHessian', true, 'logPlot', true);

%% Evaluate mismatch over the full simulation schedule 
[misfitT,~,wellSolsT] = ...
    evaluateMatch(popt,mismatchFn,trainSetup,trainPrms,trainStates,'Gradient','none');
[misfitP, ~,wellSolP,solsP] = ...
    evaluateMatch(popt,mismatchFn,predSetup,predPrms,predStates,'Gradient','none');

%% Plot well responses for training data
trainSteps = predSchedule.step.val;
fh = plotWellSols({trainWellSols,wellSolT0,wellSolsT}, trainSteps, ...
    'datasetnames', {'data','init','match'}, 'zoom', true, ...
    'field', 'qWs', 'SelectedWells', 9:10);
set(fh, 'name','Egg: GPSNet training')


%% Plot well responses for prediction case
predSteps = predSchedule.step.val;
fh = plotWellSols({predWellSols,wellSolP0, wellSolP}, predSteps, ...
    'datasetnames', {'data','init','match'}, 'zoom', true, ...
    'field', 'qWs', 'SelectedWells', 9:10);
set(fh, 'name','Egg: GPSNet prediction')

%% Plot the evolving saturation
figure
gpsNet.plotGrid(solsP{1}.s(:,1)); colormap(flipud(winter)); caxis([0 1])
for i=2:numel(solsP)
    cla
    gpsNet.plotGrid(solsP{i}.s(:,1)); drawnow; pause(0.1);
end
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
