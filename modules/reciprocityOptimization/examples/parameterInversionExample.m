%% Illustrate use Adjoint optimization for parameter inversion/reconstruction
% In this example we explain how to construct/complete a missing parameter 
% that simultanously match  water/oil production rates and bhp observation 
% from wells. The model is simple [9 9 6] grid oil water model with four wells.
%
% We construct obervations of water/oil production rates and bhp by solving the 
% model for the correct parmeter. For illustration, we propose to construct porevolume from: 
%   - water/oil production rates and bhp from each well, 
%   - Pore volume of each well cell
%   - boundary data on top and bottom faces of wells mimick seismic
%   pressure and saturation
clear all;
close  all;
mrstModule add agglom upscaling coarsegrid...
               ad-core ad-blackoil ad-props...
               optimization deckformat

time_steps     = [0.1, 0.25, 0.75, 1 ,1,  1*ones(1,20)]*day(); 
inverting_steps = 1:15;

%% Setting up the direct model
% We make a small model that consists of two different facies with
% contrasting petrophysical properties. Two injectors and two producers are
% placed in the corner at diferent depths.
gravity off;
[modelRef, W, state0] = simpleModelForInversionExample();

figure(1);
plotCellData(modelRef.G,modelRef.rock.poro,'EdgeAlpha',.5); view(3);
title('Reference model porosity')
plotWell(modelRef.G,W,'Color','k'); axis off tight
drawnow

%% Run model with exact parameter
scheduleRef = simpleSchedule(time_steps, 'W', W);
scheduleRef.step.val     = scheduleRef.step.val(1:25);
scheduleRef.step.control = scheduleRef.step.control(1:25);
Wref     = scheduleRef.control.W;
dt       = scheduleRef.step.val;
nstep    = numel(scheduleRef.step.val);
gravity off;
% For the inversion run, we create a schedule with varying controls. As
% much as this is complex scenario as long the inversion can be efficient
perturbStep = [ones(1,4), round(.5+(2:(nstep-3))/2)]';
rng(0)
scheduleRef = perturbedSimpleSchedule(dt, 'W', Wref, ...
    'pressureFac', .05, 'rateFac', .6, 'perturbStep', perturbStep);
% modelRef.toleranceCNV = 1e-6;
modelRef.useCNVConvergence = false; 
modelRef.nonlinearTolerance=1.0e-9;
modelRef.FacilityModel.toleranceWellBHP = 90000;
modelRef.FacilityModel.toleranceWellRate = 1.0e-8;
modelRef.nonlinearTolerance = 1.0e-7;
[wellSolsRef, statesRef] = simulateScheduleAD(state0, modelRef, scheduleRef);
modelRef.OutputStateFunctions{end+1} = 'ComponentTotalFlux';
figure(1)
subplot(1,2,2); cla;
plotCellData(modelRef.G, modelRef.rock.poro,'EdgeColor','none');
title('Exact PV')
plotFaces(modelRef.G, boundaryFaces(modelRef.G), 'EdgeColor', [0.4 0.4 0.4], ...
         'EdgeAlpha',.5, 'FaceColor', 'none'); view(3);
plotWell(modelRef.G, W, 'Color', 'k'); axis off tight

%%  We design simulations on the same schedule but uses different well controls. The first 
%%  modelBhpRate use the same wll controls as the standard model while modelRateBhp 
%% flip controls 
% To verify cobservation ompatibility of well output we compare the two simulators 
% for the exact parameter. 
scheduleBhpRate= scheduleRef;
modelBhpRate = modelRef;
[wellSolsBhpRate, statesBhpRate] = simulateScheduleAD(state0,modelBhpRate , scheduleBhpRate);

%% Now we create the second schedule
WW=scheduleRef.control.W;
if modelRef.G.griddim>2
    maxDepth = max(modelRef.G.faces.centroids(:,3));
    minDepth = min(modelRef.G.faces.centroids(:,3));
end

if modelRef.G.griddim<3
    maxDepth = max(modelRef.G.faces.centroids(:,2));
    minDepth = min(modelRef.G.faces.centroids(:,2));
end
[bc, bd_cells, bd_faces, well_cells] = introduceBoundaryConditionsFromSchdule(modelRef.G, statesRef, scheduleRef, maxDepth, minDepth, 'onlywells', 0);
scheduleRateBhp = flipWellControlsAndAddBC(scheduleRef, statesRef, bc);
[wellSolsRateBhp, statesRateBhp] = simulateScheduleAD(state0, modelBhpRate, scheduleRateBhp);

summary_plots2 = plotWellSols({wellSolsRef, wellSolsBhpRate} ,{scheduleRef.step.val, scheduleBhpRate.step.val},...
                            'datasetnames',{'ref. model','rate-bhp model'});
summary_plots3 = plotWellSols({wellSolsBhpRate, wellSolsRateBhp} ,{scheduleBhpRate.step.val, scheduleRateBhp.step.val},...
                            'datasetnames',{'bhp-rate model','rate-bhp model (compatibility)'});

%% Prepare parameters with corresponding bounds 
params_modelBhpRate = [];
params_modelRateBhp = [];

%% Initialize a model with an estimated porosity: here we consider max(pv) from the reference model
InitRock= modelRef.rock;
InitRock.poro = InitRock.poro.*0.0 +0.40094;
InitModel = GenericBlackOilModel(modelRef.G, InitRock, modelRef.fluid, 'gas', false);
well_data = modelRef.operators.pv(well_cells);
InitModel.operators.pv(well_cells,:)= well_data;
InitModel.toleranceCNV = 1e-8;
InitModel.nonlinearTolerance=1.0e-8;
InitModel.FacilityModel.toleranceWellBHP = 9000;
InitModel.FacilityModel.toleranceWellRate = 1.0e-8;
InitModel.FacilityModel.nonlinearTolerance =1.0e-8;
InitModel.OutputStateFunctions{end+1} = 'ComponentTotalFlux';
% We run first two simulations with compatible schedule_1 and schedule_2
% with initial parameters
modelBhpRate = InitModel;
modelRateBhp =  InitModel;
[wellSolsBhpRate, statesBhpRate] = simulateScheduleAD(state0, modelBhpRate, scheduleBhpRate);
[wellSolsRateBhp, statesRateBhp] = simulateScheduleAD(state0, modelRateBhp, scheduleRateBhp);

summary_plots3 = plotWellSols({wellSolsBhpRate, wellSolsRateBhp} ,{scheduleBhpRate.step.val, scheduleRateBhp.step.val},...
                            'datasetnames',{'bhp-rate model','rate-bhp model (initial)'});
%% parameter setting
setup_modelBhpRate = struct('model', InitModel, 'schedule', scheduleBhpRate, 'state0', state0);
setup_modelRateBhp = struct('model', InitModel, 'schedule', scheduleRateBhp, 'state0', state0);
tol = 0.01;
maxpv = min(modelRef.operators.pv)-tol;
minpv = max(modelRef.operators.pv)+tol;
params_modelBhpRate = addParameter(params_modelBhpRate,setup_modelBhpRate, 'name', 'porevolume','boxLims', [maxpv minpv]);
params_modelRateBhp = addParameter(params_modelRateBhp,setup_modelRateBhp, 'name', 'porevolume','boxLims', [maxpv minpv]);

%%  inverse model
p0 = getScaledParameterVector(setup_modelBhpRate, params_modelBhpRate);   
 weighting =  {'WaterRateWeight',  (100*stb/day)^-1, ...
 'OilRateWeight',    (100*stb/day)^-1,...
 'BHPWeight',        (95*barsa)^-1};

initBFGS = true;
itBFGS = 10;
maxIt = 30;
if initBFGS
% define objective
obj = @(model_1, states_1, schedule_1, model_2, states_2, schedule_2,states_c, doPartials, tstep, state1,state2)...
    matchObservedTwinModelsOW(model_1, states_1, schedule_1, model_2, states_2, schedule_2,states_c,bd_cells,bd_faces,well_cells,well_data, ...
    'computePartials', doPartials, 'tstep', tstep, 'state1', state1,'state2', state2, ...
    'from_states', false, weighting{:});
% provide handle to function that runs simulation for provided parameters, 
% computes objective and gradient (if requested).
objh = @(p)evaluateMatchTwinModels(p, obj, setup_modelBhpRate,setup_modelRateBhp, params_modelBhpRate, statesRef,well_cells,well_data);
% run optimization
[v, p_opt, history] = unitBoxBFGS(p0, objh, 'objChangeTol', 1e-7, 'gradTol', 1e-7, ...
                                 'maxIt', itBFGS, 'lbfgsStrategy', 'dynamic', 'lbfgsNum', 5);

%% Re-run simulation for optimal parameters for full time-horizon
setup_opt = updateSetupFromScaledParameters(setup_modelBhpRate, params_modelBhpRate, p_opt);
setup_opt.schedule.step = scheduleRef.step; % full time-horizon
[wellSols_opt, states_opt] = simulateScheduleAD(setup_opt.state0, setup_opt.model, setup_opt.schedule);
end



%% Model calibration Levenberg-Marquard (using full Jacobian)
obj2 = @(model_1, states_1, schedule_1, model_2, states_2, schedule_2,states_c, compDer, tstep, state1,state2)...
    matchObservedTwinModelsOW(model_1, states_1, schedule_1, model_2, states_2, schedule_2,states_c,bd_cells,bd_faces,well_cells,well_data, ...
    'computePartials', compDer, 'tstep', tstep, 'state1', state1,'state2', state2, ...
    'from_states', false, weighting{:}, 'mismatchSum', false);
objh2 = @(p) evaluateMatchTwinModelsSummands(p, obj2, setup_modelBhpRate,setup_modelRateBhp, params_modelBhpRate, statesRef,well_cells,well_data);
% The calibration can be improved by taking a large number of iterations,
% but here we set a limit of 30 iterations
[v2, p_opt2, history2] = unitBoxLM(p_opt, objh2, 'maxIt', maxIt);


%% Re-run simulation for optimal parameters for full time-horizon
setup_opt2 = updateSetupFromScaledParameters(setup_modelBhpRate, params_modelBhpRate, p_opt2);
setup_opt2.schedule.step = scheduleRef.step; % full time-horizon
[wellSols_opt2, states_opt2] = simulateScheduleAD(setup_opt2.state0, setup_opt2.model, setup_opt2.schedule);

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
summary_plots3 = plotWellSols({wellSolsBhpRate, wellSolsRateBhp} ,{scheduleBhpRate.step.val, scheduleRateBhp.step.val},...
                            'datasetnames',{'fine scale model','initial upscaled model'});