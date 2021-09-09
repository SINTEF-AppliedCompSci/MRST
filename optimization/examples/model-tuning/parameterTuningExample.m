%% Illustrate use Adjoint optimization for parameter tuning
% In this example we explain how to tune parameters of a coarse model to 
% best match the output of a fine model. 
%
% We tune the coarse model to match water/oil production rates and bhp. For
% illustration, we indroduce and tune many parameters: 
%   -Transmisibility of each internal faces, 
%   -Pore volume of each cell
%   -Well index for each well 
%   -Six relperm scalers per cell

mrstModule add agglom upscaling coarsegrid...
               ad-core ad-blackoil ad-props...
               optimization deckformat

time_steps     = [1 ,1 ,2 5 10 , 10*ones(1,10)]*day(); 
training_steps = 1:10;
partition      = [3,3,2];

%% Setting up the fine-scale model
% We make a small model that consists of two different facies with
% contrasting petrophysical properties. Two injectors and two producers are
% placed in the corner at diferent depths.
[model_f, W, state0] = simpleModelForTuningExample();

figure(1), clf, subplot(1,2,1)
plotCellData(model_f.G,model_f.rock.poro,'EdgeAlpha',.5); view(3);
title('Fine-scale model porosity')
plotWell(model_f.G,W,'Color','k'); axis off tight
drawnow

%% Run simulation for fine-scale model
schedule = simpleSchedule(time_steps, 'W', W);
problem  = packSimulationProblem(state0, model_f, schedule, 'model_fine_scale');
[ok, status] = simulatePackedProblem(problem);
[wellSols_f, states_f] = getPackedSimulatorOutput(problem);

%% Coarse-scale model
% We make a coarse grid defined by partition, and perform a simple
% upscaling to obtain a coarse model
p  = partitionCartGrid(model_f.G.cartDims, partition );
model_c = upscaleModelTPFA(model_f, p);
% We want to include rel-perm scalers as tunabale parameters, so include
% these for the coarse model. Scalers have no effect for the initial coarse 
% model (they are set equal to the ones given by the rel-perm curves). 
pts = model_c.fluid.krPts;
scaling = {'SWL',   pts.w(1), 'SWCR', pts.w(2), 'SWU', pts.w(3), ...
           'SOWCR', pts.o(2), 'KRW',  pts.w(4), 'KRO', pts.o(4)};
model_c = imposeRelpermScaling(model_c, scaling{:});
% use a tighter tollerance for improved gradient accuracy
model_c.toleranceCNV = 1e-6;
% perform a simple upscaling of the schedule for the training runs
schedule_training   = simpleSchedule(time_steps(training_steps), 'W', W);
schedule_training_c = upscaleSchedule(model_c, schedule_training);


figure(1)
subplot(1,2,2); cla;
plotCellData(model_c.G, model_c.rock.poro,'EdgeColor','none');
title('Coarse-scale model porosity')
plotFaces(model_c.G, boundaryFaces(model_c.G), 'EdgeColor', [0.4 0.4 0.4], ...
         'EdgeAlpha',.5, 'FaceColor', 'none'); view(3);
plotWell(model_f.G, W, 'Color', 'k'); axis off tight

 %% Simulate initial upscaled coarse model for full time
state0_c   = upscaleState(model_c, model_f, state0);
schedule_c = upscaleSchedule(model_c, schedule);
[wellSols_c, states_c] = simulateScheduleAD(state0_c, model_c, schedule_c);
summary_plots = plotWellSols({wellSols_f, wellSols_c} ,{schedule.step.val, schedule_c.step.val},...
                            'datasetnames',{'fine scale model','initial upscaled model'});
drawnow

%% Prepare parameters with corresponding bounds 
setup = struct('model', model_c, 'schedule', schedule_training_c, 'state0', state0_c);
params = [];
% Well, porevolume and transmisibility
params = addParameter(params,setup, 'name', 'conntrans',       'relativeLimits', [.25 4], 'scaling', 'log');
params = addParameter(params,setup, 'name', 'porevolume',      'relativeLimits', [.5  2]);
params = addParameter(params,setup, 'name', 'transmissibility','relativeLimits', [.25 4], 'scaling', 'log');
% Rel-perm scalers
params = addParameter(params,setup, 'name', 'swl',  'boxLims',[0.0 0.2]);
params = addParameter(params,setup, 'name', 'swcr', 'boxLims',[0.0 0.2]);
params = addParameter(params,setup, 'name', 'sowcr','boxLims',[0.0 0.2]);
params = addParameter(params,setup, 'name', 'swu',  'boxLims',[0.8 1.0]);
params = addParameter(params,setup, 'name', 'kro',  'boxLims',[0.8 1]);
params = addParameter(params,setup, 'name', 'krw',  'boxLims',[0.8 1]);


%% Coarse model tuning
p0 = getScaledParameterVector(setup, params);
% Weight production rates/bhp relative to their magnitudes
weighting =  {'WaterRateWeight',  (5/day)^-1, ...
              'OilRateWeight',    (5/day)^-1,...
              'BHPWeight',        (5*barsa)^-1};            
% define objective
obj = @(model, states, schedule, states_ref, doPartials, tstep, state)...
    matchObservedOW(model, states, schedule, states_ref, ...
    'computePartials', doPartials, 'tstep', tstep, 'state', state, ...
    'from_states', false, weighting{:});
% provide handle to function that runs simulation for provided parameters, 
% computes objective and gradient (if requested).
objh = @(p)evaluateMatch(p, obj, setup, params, states_f);
% run optimization
[v, p_opt, history] = unitBoxBFGS(p0, objh, 'objChangeTol', 1e-7, 'gradTol', 1e-4, ...
                                 'maxIt', 25, 'lbfgsStrategy', 'dynamic', 'lbfgsNum', 10);

%% Re-run simulation for optimal parameters for full time-horizon
setup_opt = updateSetupFromScaledParameters(setup, params, p_opt);
setup_opt.schedule.step = schedule.step; % full time-horizon
[wellSols_opt, states_opt] = simulateScheduleAD(setup_opt.state0, setup_opt.model, setup_opt.schedule);
                                                        
try
    figure(summary_plots.Number)
catch
    summary_plots = figure;
end
ts = schedule.step.val;
plotWellSols({wellSols_f, wellSols_c, wellSols_opt, wellSols_f(1:10)},...
             {ts,         ts,         ts,           ts(1:10)},...
              'datasetnames',{'reference','initial coarse', 'tuned coarse', 'training points'},...
              'linestyles',{'-', '-', '-', 'o'},...
              'figure',summary_plots.Number)

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
