%% Norne: calibration of a network model formed from a coarse grid
% This example demonstrates how to set up a CGNet model with 3D type
% graph topology and corresponding parameters initialized based on a coarse
% partition of the full corner-point grid of the Norne field model.
%
% To calibrate the model, we use the BroydenFletcherGoldfarbShanno
% (BFGS) algorithm. This is an iterative line-search method that gradually
% improves an approximation to the Hessian matrix of the mismatch function,
% obtained only from adjoint gradients via a generalized secant method.
%
% This example was first introduced in MRST 2021b.
mrstModule add ad-core ad-blackoil deckformat ...
               agglom upscaling coarsegrid book ...
               mrst-gui ad-props incomp optimization...
               network-models test-suite linearsolvers 

%% Setup 3D reference model
% The reference model is a single stochastic realization of the Norne field
% model taken from the example-suit module. Compared with the real field,
% this simulation case has simpler fluid description (an oil-water model)
% and an idealized field development plan consisting of a simple pattern of
% eleven vertical wells that run under constant bhp or rate controls.
%
% We construct two different schedules: a plausible schedule used to
% produce the field and a schedule with random perturbations in the well
% controls to train the data-driven model

% True schedule, which we seek to reproduce
predCase  = TestCase('norne_simple_wo');
predProbl = predCase.getPackedSimulationProblem();
%clearPackedSimulatorOutput(trueProb)
simulatePackedProblem(predProbl);

[predWellSols, predStates] = getPackedSimulatorOutput(predProbl);
predModel     = predCase.model;
predSchedule  = predProbl.SimulatorSetup.schedule;
Wtrue         = predSchedule.control.W;

% Random schedule
trainCase  = makeRandomTraining(predCase, 0.25, 0.05, false);
trainProbl = trainCase.getPackedSimulationProblem();

%clearPackedSimulatorOutput(trainProb)
simulatePackedProblem(trainProbl);

[trainWellSols, trainStates] = getPackedSimulatorOutput(trainProbl);
trainModel     = trainCase.model;
trainSchedule  = trainProbl.SimulatorSetup.schedule;
Wtrain         = trainSchedule.control.W;

plotWellSols({predWellSols, trainWellSols}, ...
    {trainSchedule.step.val, predSchedule.step.val},...
    'datasetnames',{'reference','training'}, ...
    'zoom', true, 'field', 'qOs', 'SelectedWells', 7);

%% Coarse-scale model and initial state
% We make a coarse grid defined by a uniform 6 x 8 x 1 partition and
% perform a simple upscaling to obtain a coarse model
q = processPartition(trainModel.G, partitionUI(trainModel.G,[5 6 1]));
q = compressPartition(q);
cModel = upscaleModelTPFA(trainModel, q,'transFromRock',false);

% We want to include rel-perm scaling as tunabale parameters, so include
% these for the coarse model. These parameters have no effect for the
% initial coarse model (they are set equal to the ones given by the
% rel-perm curves).
pts = cModel.fluid.krPts;
scaling = {'SWL',   pts.w(1,1), 'SWCR', pts.w(1,2), 'SWU', pts.w(1,3), ...
           'SOWCR', pts.o(1,2), 'KRW',  pts.w(1,4), 'KRO', pts.o(1,4)};
cModel = imposeRelpermScaling(cModel, scaling{:});
cModel = cModel.setupOperators();

cModel.toleranceCNV = 1e-6;  % tighter tolerance to improve gradient accuracy

cState0 = upscaleState(cModel, trainModel, trainCase.state0);

%% Specify training schedule and parameters to be matched
% We use the random schedule for training. In this setup, we use all pore
% volumes, transmissibilities, and well connections in the coarse grid as
% calibration parameters.
cTrainSched = upscaleSchedule(cModel, trainSchedule);
cTrainProbl = struct('model', cModel, 'schedule', cTrainSched, 'state0', cState0);
cPredSched  = upscaleSchedule(cModel, predSchedule);
cPredProbl  = struct('model', cModel, 'schedule', cPredSched, 'state0', cState0);

config = {
    ...%name      include     scaling    boxlims  lumping   subset  relativeLimits
    'porevolume',       1,   'linear',       [],    [],      [],    [.001 4]
    'conntrans',        1,   'log',          [],    [],      [],    [.001 100]
    'transmissibility', 1,   'log'   ,       [],    [],      [],    [.001 100]};
trainPrms = []; predPrms = [];
for k = 1:size(config,1)
    if config{k, 2} == 0, continue, end
    trainPrms = addParameter(trainPrms, cTrainProbl, ...
        'name',    config{k,1}, 'scaling', config{k,3}, ...
        'boxLims', config{k,4}, 'lumping', config{k,5}, ...
        'subset',  config{k,6}, 'relativeLimits',config{k,7});
    predPrms = addParameter(predPrms, cPredProbl, ...
        'name',    config{k,1}, 'scaling', config{k,3}, ...
        'boxLims', config{k,4}, 'lumping', config{k,5}, ...
        'subset',  config{k,6}, 'relativeLimits',config{k,7});
end

%% Plot the CGNet model
fig1 = figure;
subplot(2,2,1)
G = predCase.getVisualizationGrid();
plotCellData(G, predCase.model.rock.poro,'EdgeColor','none');
plotWell(G, predCase.schedule.control(1).W,'color','k','FontSize',10);
view(85,65); 
axis tight off; set(gca,'Clipping',false); zoom(1.2);

subplot(2,2,3), ax=gca;
nq = max(q);
G = predCase.getVisualizationGrid();
colormap(gca,tatarizeMap(nq));
explosionView(G,q);
set(ax.Children,'EdgeAlpha',.1); view(85,75); 
axis tight off; set(gca,'Clipping',false); zoom(1.2);

subplot(2,2,4), ax=gca;
CG = cModel.G;
N = getNeighbourship(CG);
A = getConnectivityMatrix(N,true,CG.cells.num);
network = graph(A-eye(CG.cells.num));
pg = plot(network,'LineWidth',2, ...
    'XData', CG.cells.centroids(:,1), ...
    'YData', CG.cells.centroids(:,2), ...
    'ZData', CG.cells.centroids(:,3));
labelnode(pg,1:nq,'');
cells=vertcat(cPredSched.control(1).W.cells);
names=rldecode({cPredSched.control(1).W.name},...
    vertcat(arrayfun(@(x) numel(x.cells), cPredSched.control(1).W)),2);
labelnode(pg,cells,names);
set(ax.Children,'NodeFontSize',10,'NodeFontWeight','bold');
plotGrid(G,'FaceColor','none','EdgeAlpha',.05);
view(85,75); 
axis tight off; set(gca,'Clipping',false); zoom(1.2);

subplot(2,2,2), ax=gca;
pg = plot(network,'LineWidth',0.5,'Layout','circle');
labelnode(pg,1:nq,'');
labelnode(pg,cells,names);
set(ax.Children,'NodeFontSize',10,'NodeFontWeight','bold');
axis equal tight off; set(gca,'Clipping',false); zoom(1.2);

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
                   'state', state, 'from_states', false);

%% Model calibration
% Calibrate the model using the BFGS method. This is a computationally
% expensive operation that may run for several hours if you choose a large
% number of iterations. Here, we therefore only apply 10 iterations. To get
% a good match, you should increase this number to 100+
pinit = getScaledParameterVector(cTrainProbl, trainPrms);
objh = @(p) evaluateMatch(p,mismatchFn,cTrainProbl,trainPrms,trainStates);
[v, popt, history] = unitBoxBFGS(pinit, objh, 'objChangeTol', 1e-8, ...
    'gradTol', 1e-5, 'maxIt',10, 'lbfgsStrategy', 'dynamic', ...
    'lbfgsNum', 5, 'outputHessian',true, 'logPlot', true);

%% Evaluate mismatch over the full simulation schedule 
[misfitP,~,wellSolP] = ...
    evaluateMatch(popt,mismatchFn,cPredProbl,predPrms,predStates,'Gradient','none');
[misfitT,~,wellSolT] = ...
    evaluateMatch(popt,mismatchFn,cTrainProbl,trainPrms,trainStates,'Gradient','none');

%% Plot well curves
fh = plotWellSols({trainWellSols,wellSolT}, ...
    {trainSchedule.step.val,cTrainSched.step.val}, ...
    'datasetnames', {'train','match'}, 'zoom', true, ...
    'field', 'qWs', 'SelectedWells', 7:11);
set(fh, 'name','Norne')

%%
fh = plotWellSols({predWellSols,wellSolP}, ...
    {predSchedule.step.val,cPredSched.step.val}, ...
    'datasetnames', {'reference','predicted'}, 'zoom', true, ...
    'field', 'qOs', 'SelectedWells', 7);
set(fh, 'name','Norne')

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
                                     