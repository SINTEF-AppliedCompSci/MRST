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

clear all;
close all;
mrstModule add agglom upscaling coarsegrid...
               ad-core ad-blackoil ad-props...
               optimization deckformat
mrstModule add ad-core ad-props ad-blackoil spe10
time_steps     = [0.1, 0.1, 0.25, 0.5,1 ,2 5 10 , 10*ones(1,10)]*day(); 
training_steps = 1:18;
partition      = [3,3,2];

%% Setting up the fine-scale model
% We make a small model that consists of two different facies with
% contrasting petrophysical properties. Two injectors and two producers are
% placed in the corner at diferent depths.
[model_f, W, state0] = simpleModelForTuningExample();
model_f.OutputStateFunctions{end + 1} = 'Mobility';
model_f.nonlinearTolerance=1.0e-9;
model_f.toleranceCNV = 1.0e-10;
model_f.toleranceMB = 1.0e-7;
model_f.FacilityModel.toleranceWellBHP = 1;
model_f.FacilityModel.toleranceWellRate = 1.0e-11;
model_f.FacilityModel.nonlinearTolerance =1.0e-9;
figure(1), clf, subplot(1,2,1)
plotCellData(model_f.G,model_f.rock.poro,'EdgeAlpha',.5); view(3);
title('Fine-scale model porosity')
plotWell(model_f.G,W,'Color','k'); axis off tight
drawnow

%% Run simulation for fine-scale model
scheduleWellNoBc = simpleSchedule(time_steps, 'W', W);
problemWellNoBc  = packSimulationProblem(state0, model_f, scheduleWellNoBc, 'model_fine_scale_Well');
problemWellNoBc.Modules(end-1:end) = [];
simulatePackedProblem(problemWellNoBc);%, 'restartStep',1);
[wellSols_f_WellNoBc, states_f_WellNoBc] = getPackedSimulatorOutput(problemWellNoBc);

%% Run another fine-scale with Dirichlet boundary conditions
%scheduleWellBc = scheduleWellNoBc;
bdf_f = boundaryFaces(model_f.G);
cells_bdf= sum(model_f.G.faces.neighbors(bdf_f,:),2);
cellNo = rldecode(1:model_f.G.cells.num, diff(model_f.G.cells.facePos), 2) .';
for i = 1:length(time_steps)
%     scheduleWellBc.control(i).W = scheduleWellNoBc.control.W;
    fpress     =  states_f_WellNoBc{i}.pressure(cells_bdf);
   % scheduleWellBc.step.control(i) = i;
    bc{i} =  addBC([], bdf_f, 'pressure',        ...
            fpress,'sat', [states_f_WellNoBc{i}.s(cells_bdf,:)]);
end
scheduleWellBc = flipWellControlsAndAddBC(scheduleWellNoBc, states_f_WellNoBc, []);

problemWellBc  = packSimulationProblem(state0, model_f, scheduleWellBc, 'model_fine_scale_WellBc');
simulatePackedProblem(problemWellBc);%,'restartStep',1);
[wellSols_f_WellBc, states_f_WellBc] = getPackedSimulatorOutput(problemWellBc);
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
%model_c.toleranceCNV = 1e-6;
% perform a simple upscaling of the schedule for the training runs
schedule_training_WellNoBc   = simpleSchedule(time_steps(training_steps), 'W', W);
schedule_training_c_WellNoBc = upscaleSchedule(model_c, schedule_training_WellNoBc,   'bcUpscaleMethod', 'idw', 'wellUpscaleMethod', 'recompute');
% another schedule
schedule_training_WellBc.control = scheduleWellBc.control;
for i = training_steps
schedule_training_WellBc.step.val(i) = scheduleWellBc.step.val(i);
schedule_training_WellBc.step.control(i) = scheduleWellBc.step.control(i);
end
schedule_training_c_WellBc = upscaleSchedule(model_c, schedule_training_WellBc, 'bcUpscaleMethod', 'idw', 'wellUpscaleMethod', 'recompute');



figure(1)
subplot(1,2,2); cla;
plotCellData(model_c.G, model_c.rock.poro,'EdgeColor','none');
title('Coarse-scale model porosity')
plotFaces(model_c.G, boundaryFaces(model_c.G), 'EdgeColor', [0.4 0.4 0.4], ...
         'EdgeAlpha',.5, 'FaceColor', 'none'); view(3);
plotWell(model_f.G, W, 'Color', 'k'); axis off tight

 %% Simulate initial upscaled coarse model for full time
state0_c   = upscaleState(model_c, model_f, state0);
% state0_c.wellSol(1).pressure =states_f_WellNoBc{1}.wellSol(1).bhp;
% state0_c.wellSol(2).pressure =states_f_WellNoBc{1}.wellSol(2).bhp;
% state0_c.wellSol(3).pressure =states_f_WellNoBc{1}.wellSol(3).bhp;
% state0_c.wellSol(4).pressure =states_f_WellNoBc{1}.wellSol(4).bhp;
schedule_c_WellNoBc = upscaleSchedule(model_c, scheduleWellNoBc, 'bcUpscaleMethod', 'idw', 'wellUpscaleMethod', 'recompute');
schedule_c_WellBc = upscaleSchedule(model_c, scheduleWellBc,   'bcUpscaleMethod', 'idw', 'wellUpscaleMethod', 'recompute');
[wellSols_c_WellNoBc, states_c_WellNoBc] = simulateScheduleAD(state0_c, model_c, schedule_c_WellNoBc);
[wellSols_c_WellBc, states_c_WellBc] = simulateScheduleAD(state0_c, model_c, schedule_c_WellBc);

summary_plots = plotWellSols({wellSols_f_WellNoBc, wellSols_c_WellNoBc, wellSols_f_WellBc, wellSols_c_WellBc} ,{scheduleWellNoBc.step.val, schedule_c_WellNoBc.step.val, scheduleWellBc.step.val, schedule_c_WellBc.step.val},...
                            'datasetnames',{'fine scale WellCtr model','initial upscaled WellCtr model', 'fine scale WellBcCtr model','initial upscaled WellBcCtr model'});
drawnow

%% Prepare parameters with corresponding bounds 
setup_WellNoBc = struct('model', model_c, 'schedule', schedule_training_c_WellNoBc, 'state0', state0_c);
setup_WellBc = struct('model', model_c, 'schedule', schedule_training_c_WellBc, 'state0', state0_c);

params_WellNoBc = [];
params_WellBc = [];
% Well, porevolume and transmisibility
params_WellNoBc = addParameter(params_WellNoBc,setup_WellNoBc, 'name', 'conntrans',       'relativeLimits', [.25 4], 'scaling', 'log');
params_WellNoBc = addParameter(params_WellNoBc,setup_WellNoBc, 'name', 'porevolume',      'relativeLimits', [.5  2]);
params_WellNoBc = addParameter(params_WellNoBc,setup_WellNoBc, 'name', 'transmissibility','relativeLimits', [.25 4], 'scaling', 'log');
% Rel-perm scalers
params_WellNoBc = addParameter(params_WellNoBc,setup_WellNoBc, 'name', 'swl',  'boxLims',[0.0 0.2]);
params_WellNoBc = addParameter(params_WellNoBc,setup_WellNoBc, 'name', 'swcr', 'boxLims',[0.0 0.2]);
params_WellNoBc = addParameter(params_WellNoBc,setup_WellNoBc, 'name', 'sowcr','boxLims',[0.0 0.2]);
params_WellNoBc = addParameter(params_WellNoBc,setup_WellNoBc, 'name', 'swu',  'boxLims',[0.8 1.0]);
params_WellNoBc = addParameter(params_WellNoBc,setup_WellNoBc, 'name', 'kro',  'boxLims',[0.8 1]);
params_WellNoBc = addParameter(params_WellNoBc,setup_WellNoBc, 'name', 'krw',  'boxLims',[0.8 1]);


% Well, porevolume and transmisibility
params_WellBc = addParameter(params_WellBc,setup_WellBc, 'name', 'conntrans',       'relativeLimits', [.25 4], 'scaling', 'log');
params_WellBc = addParameter(params_WellBc,setup_WellBc, 'name', 'porevolume',      'relativeLimits', [.5  2]);
params_WellBc = addParameter(params_WellBc,setup_WellBc, 'name', 'transmissibility','relativeLimits', [.25 4], 'scaling', 'log');
% Rel-perm scalers
params_WellBc = addParameter(params_WellBc,setup_WellBc, 'name', 'swl',  'boxLims',[0.0 0.2]);
params_WellBc = addParameter(params_WellBc,setup_WellBc, 'name', 'swcr', 'boxLims',[0.0 0.2]);
params_WellBc = addParameter(params_WellBc,setup_WellBc, 'name', 'sowcr','boxLims',[0.0 0.2]);
params_WellBc = addParameter(params_WellBc,setup_WellBc, 'name', 'swu',  'boxLims',[0.8 1.0]);
params_WellBc = addParameter(params_WellBc,setup_WellBc, 'name', 'kro',  'boxLims',[0.8 1]);
params_WellBc = addParameter(params_WellBc,setup_WellBc, 'name', 'krw',  'boxLims',[0.8 1]);


%% Coarse model tuning
p0 = getScaledParameterVector(setup_WellNoBc, params_WellNoBc);
% Weight production rates/bhp relative to their magnitudes
weighting =  {'WaterRateWeight',  [], ...
              'OilRateWeight',    [],...
              'BHPWeight',        []};


obj = @(model_1, states_1, schedule_1, model_2, states_2, schedule_2,states_c, doPartials, tstep, state1,state2)...
    matchObservedTwinModelsOW(model_1, states_1, schedule_1, model_2, states_2, schedule_2,states_c,[],[],[],[], ...
    'computePartials', doPartials, 'tstep', tstep, 'state1', state1,'state2', state2, ...
    'from_states', false, weighting{:});

% % define objective
% obj = @(model, states, schedule, states_ref, doPartials, tstep, state)...
%     matchObservedOW(model, states, schedule, states_ref, ...
%     'computePartials', doPartials, 'tstep', tstep, 'state', state, ...
%     'from_states', false, weighting{:});
% provide handle to function that runs simulation for provided parameters, 
% computes objective and gradient (if requested).
% objh = @(p)evaluateMatch(p, obj, setup_WellNoBc, params_WellNoBc, states_f_WellNoBc);

objh = @(p)evaluateMatchTwinModels(p, obj, setup_WellNoBc,setup_WellBc, params_WellNoBc, states_f_WellBc,[],[]);
% run optimization
[v, p_opt, history] = unitBoxBFGS(p0, objh, 'objChangeTol', 1e-7, 'gradTol', 1e-4, ...
                                 'maxIt', 25, 'lbfgsStrategy', 'dynamic', 'lbfgsNum', 10);

%% Re-run simulation for optimal parameters for full time-horizon
setup_opt_WellNoBc = updateSetupFromScaledParameters(setup_WellNoBc, params_WellNoBc, p_opt);
setup_opt_WellNoBc.schedule.step = scheduleWellNoBc.step; % full time-horizon
[wellSols_opt_WellNoBc, states_opt_WellNoBc] = simulateScheduleAD(setup_opt_WellNoBc.state0, setup_opt_WellNoBc.model, setup_opt_WellNoBc.schedule);
  
%% Re-run simulation for optimal parameters for full time-horizon
setup_opt_WellBc = updateSetupFromScaledParameters(setup_WellBc, params_WellBc, p_opt);
setup_opt_WellBc.schedule.step = scheduleWellBc.step; % full time-horizon
[wellSols_opt_WellBc, states_opt_WellBc] = simulateScheduleAD(setup_opt_WellBc.state0, setup_opt_WellBc.model, setup_opt_WellBc.schedule);

try
    figure(summary_plots.Number)
catch
    summary_plots = figure;
end
ts = scheduleWellNoBc.step.val;
plotWellSols({wellSols_f_WellNoBc, wellSols_opt_WellNoBc, wellSols_f_WellBc, wellSols_opt_WellBc},...
             {ts,         ts,         ts,           ts},...
              'datasetnames',{'reference','tuned WellNoBc', 'reference', 'tuned WellBc'},...
              'linestyles',{'-', '-.', '-', '-.'},...
              'figure',summary_plots.Number)

%% Copyright Notice
%
% <html>
% <p><font size="-1">
% Copyright 2009-2024 SINTEF Digital, Mathematics & Cybernetics.
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
