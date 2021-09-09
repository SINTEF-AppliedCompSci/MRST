%% Norne: parameter tuning of a very coarse upscaled model 
% This example demonstrates how to set up and initialize a model based on 
% a coarse partition of the full corner-point grid of the Norne field model.
%
% To calibrate the coarse model parameters, we optimize a mismatch-function 
% where the optimization utilize parameter sensitivities obtained by adjoint
% simulations.
mrstModule add ad-core ad-blackoil deckformat ...
               agglom upscaling coarsegrid...
               mrst-gui ad-props incomp optimization...
               example-suite linearsolvers 

%% Setup reference model
% The reference model is a single stochastic realization of the Norne field
% model taken from the example-suit module. Compared with the real field,
% this simulation case has simpler fluid description (an oil-water model)
% and an idealized field development plan consisting of a simple pattern of
% eleven vertical wells that run under constant bhp or rate controls.

example = MRSTExample('norne_simple_wo');

% % Modify the setup to use longer time steps in order to speed up script
% execution
totTime  = sum(example.schedule.step.val);
nstep    = 48;
example.schedule.step = struct('val',     ones(nstep,1)*(totTime/nstep), ...
                               'control', ones(nstep, 1));
% ix  = repmat(1:numel(schedule.step.val), [4 1]);
% schedule.step.val     = schedule.step.val(ix(:))/4;
% schedule.step.control = schedule.step.control(ix(:));
% problem.SimulatorSetup.schedule = schedule;
% example.schedule = schedule;
problem  = example.getPackedSimulationProblem();
% Simulate
simulatePackedProblem(problem);

[wellSolsRef, statesRef] = getPackedSimulatorOutput(problem);
modelRef    = problem.SimulatorSetup.model;
scheduleRef = problem.SimulatorSetup.schedule;
Wref        = scheduleRef.control.W;

% Plot
example.plot(statesRef, 'step_index', numel(statesRef))

%% Coarse-scale model
% We make a coarse grid defined by a uniform 6 x 8 x 1 partition 
blockIx = partitionUI(modelRef.G, [6, 8, 1]);
% The fourth layer of the fine model grid is mosly inactive, so several of
% the resulting coarse blocks are not connected. We process the partition
% futher to split disconnected blocks into new blocks
blockIx = processPartition(modelRef.G, blockIx);
blockIx = compressPartition(blockIx);
% Perform a simple upscaling to obtain a coarse model
cModel = upscaleModelTPFA(modelRef, blockIx);

% We want to include rel-perm scaling as tunabale parameters, so include
% these for the coarse model. These parameters have no effect for the
% initial coarse model (they are set equal to the ones given by the
% rel-perm curves).
pts = cModel.fluid.krPts;
scaling = {'SWL',   pts.w(1), 'SWCR', pts.w(2), 'SWU', pts.w(3), ...
           'SOWCR', pts.o(2), 'KRW',  pts.w(4), 'KRO', pts.o(4)};
cModel = imposeRelpermScaling(cModel, scaling{:});
cModel.toleranceCNV = 1e-6;  % tighter tolerance to improve gradient accuracy

%% Plot comparison of simulation oputput from the coarse and fine model
figure('position',[100 100 1000 400])
axes('position',[.02 .05 .48 .9]);
plotCellData(modelRef.G, modelRef.rock.poro, 'EdgeAlpha',.1); 
view(174,60);
title('Fine-scale model porosity')
plotWell(modelRef.G, Wref, 'Color', 'k', 'FontSize', 10); axis off tight
mrstColorbar(modelRef.rock.poro,'South'); cax = caxis;

axes('position',[.5 .05 .48 .9]);
plotCellData(cModel.G, cModel.rock.poro, 'EdgeColor', 'none');
title('Coarse-scale model porosity')
plotFaces(cModel.G, boundaryFaces(cModel.G), 'EdgeColor', [0.4 0.4 0.4], ...
         'EdgeAlpha',.5, 'FaceColor', 'none');
view(174,60); caxis(cax);
plotWell(modelRef.G, Wref, 'Color', 'k', 'FontSize', 10); axis off tight
mrstColorbar(cModel.rock.poro,'South');

%% Simulate initial upscaled coarse model for full time
cState0   = upscaleState(cModel, modelRef, example.state0);
cSchedule = upscaleSchedule(cModel, scheduleRef);
[cWellSols, cStates] = simulateScheduleAD(cState0, cModel, cSchedule);
plotWellSols({wellSolsRef, cWellSols}, ...
    {scheduleRef.step.val, cSchedule.step.val},...
    'datasetnames',{'fine scale model','initial upscaled model'}, ...
    'zoom', true, 'field', 'qOs', 'SelectedWells', 7);

%% Specify training schedule and parameters to be matched
% We use the first half of the given data for training. In this setup, we
% use all pore volumes, transmissibilities, and well connections in the
% coarse grid as calibration parameters.
trainSteps = 1:round(numel(scheduleRef.step.val)/2);
timeSteps  = scheduleRef.step.val(trainSteps); 
trainSched = upscaleSchedule(cModel, simpleSchedule(timeSteps, 'W', Wref));
setup_init = struct('model', cModel, 'schedule', trainSched, 'state0', cState0);

pv = cModel.operators.pv;
% set up 'matrix' for parameter options for easier editing. The specific
% limits set for the various parameters, influences the tuning/optimization 
% procedure to a large extent
config = {
    ...%name      include     scaling     boxlims            lumping   subset  relativeLimits  
    'porevolume',       1,   'linear',    [.01*pv, 1.5*pv],     [],     [],      []   
    'conntrans',        1,   'log',       [],                   [],     [],      [1e-2, 1e2]      
    'transmissibility', 1,   'log'   ,    [],                   [],     [],      [1e-2  1e2]  
    'swl',              1,   'linear',    [0 .3],               [],     [],      []
    'swcr',             1,   'linear',    [0 .4],               [],     [],      []
    'swu',              1,   'linear',    [.7 1],               [],     [],      []
    'sowcr',            1,   'linear',    [0 .4],               [],     [],      []
    'krw',              1,   'linear',    [.5 1.5],             [],     [],      []
    'kro',              1,   'linear',    [.5 1.5],             [],     [],      []
    'sw',               0,   'linear',    [0 .4 ],              [],     [],      []
    'pressure'          0,   'linear',    [],                   [],     [],      [.8 1.2]};
parameters = [];
for k = 1:size(config,1)
    if config{k, 2} == 0, continue, end
    parameters = addParameter(parameters, setup_init, ...
        'name',    config{k,1}, 'scaling', config{k,3}, ...
        'boxLims', config{k,4}, 'lumping', config{k,5}, ...
        'subset',  config{k,6}, 'relativeLimits',config{k,7});
end

%% Define the mismatch function
% The mismatch function is defined as a function handle to a library
% function from the optimization module that computes the mismatch between
% a given simulation and a reference state. For an oil-water system, the
% match is computed based on three quantities (water/oil rate and bhp) and
% these should be given an associated weight. Weights on the order the
% reciprocal of rate magnitudes/pressure variations result in a properly 
% scaled mismatch function 
weighting  = {'WaterRateWeight',  1/(1e4/day), ...
              'OilRateWeight',    1/(1e4/day), ...
              'BHPWeight',        1/(500*barsa)};
mismatchFn = @(model, states, schedule, states_ref, tt, tstep, state) ...
    matchObservedOW(model, states, schedule, states_ref,...
                   'computePartials', tt, 'tstep', tstep, weighting{:},...
                   'state', state, 'from_states', false);

%% Model calibration
pvec = getScaledParameterVector(setup_init, parameters);
objh = @(p) evaluateMatch(p, mismatchFn, setup_init, parameters, statesRef);
% The calibration can be improved by taking a large number of iterations,
% but here we set a limit of 30 iterations
[v, p_opt, history] = unitBoxBFGS(pvec, objh, 'objChangeTol', 1e-4, ...
    'maxIt', 30, 'lbfgsStrategy', 'dynamic', 'lbfgsNum', 5);

%% Create a new coarse model setup with the optimal parameters, 
%% and rerun the simulation for full time horizon
setup_opt = updateSetupFromScaledParameters(setup_init, parameters, p_opt); 
% reset time-steps to full schedule
setup_opt.schedule.step = scheduleRef.step;
[wellSols_opt, states_opt] = simulateScheduleAD(setup_opt.state0, setup_opt.model, setup_opt.schedule);

%% Plot well curves for reference fine, initial coarse and tuned coarse models 
fh = plotWellSols({wellSolsRef,cWellSols,wellSols_opt}, ...
    repmat({scheduleRef.step.val}, 1, 3), ...
    'datasetnames', {'reference','initial','optimized'}, 'zoom', true, ...
    'field', 'qOs', 'SelectedWells', 7);

%% Plot the pore volume updates
% fetch pore volume differences in initial and tuned coarse models
dpv = setup_opt.model.operators.pv - setup_init.model.operators.pv;
figure
plotCellData(cModel.G, dpv, 'EdgeColor','none');
plotFaces(cModel.G, boundaryFaces(cModel.G), 'EdgeColor', [0.4 0.4 0.4], ...
         'EdgeAlpha',.5, 'FaceColor', 'none');
view(174,60);
plotWell(modelRef.G, Wref, 'Color', 'k', 'FontSize', 10); axis off tight
colorbar('south');

%% Copyright Notice
%
% <html>
% <p><font size="-1">
% Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.
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
