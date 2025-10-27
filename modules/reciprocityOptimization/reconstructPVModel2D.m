clear all;
mrstModule add ad-core ad-blackoil ad-props optimization spe10
setupModel2D
gravity off;
% Create model-object of class TwoPhaseOilWaterModel
modelRef = GenericBlackOilModel(G, rock, fluid, 'gas', false);
modelRef.OutputStateFunctions(end+1) = {'ComponentTotalFlux'};
WRef     = schedule.control.W;
dt       = schedule.step.val;
nstep    = numel(schedule.step.val);

% For the training run, we create a schedule with varying controls
perturbStep = [ones(1,4), round(.5+(2:(nstep-3))/2)]';
rng(0)
scheduleRef = perturbedSimpleSchedule(dt, 'W', WRef, ...
    'pressureFac', .01, 'rateFac', .4, 'perturbStep', perturbStep);
modelRef.toleranceCNV = 1e-9;
modelRef.nonlinearTolerance=1.0e-8;
modelRef.FacilityModel.toleranceWellBHP = 1;
modelRef.FacilityModel.toleranceWellRate = 1.0e-9;
modelRef.FacilityModel.nonlinearTolerance =1.0e-9;
% Set initial state and run simulation:
state0 = initResSol(G, pRef, [0, 1]);
[wsRef, statesRef] = simulateScheduleAD(state0, modelRef, scheduleRef);


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
[wsBhpRate, statesBhpRate] = simulateScheduleAD(state0,modelBhpRate , scheduleBhpRate);

%% Now we create the second schedule
WW=scheduleRef.control.W;
if modelRef.G.griddim>2&&all(G.cells.centroids(:, 3)~= G.cells.centroids(1, 3))
    maxDepth = max(modelRef.G.faces.centroids(:,3));
    minDepth = min(modelRef.G.faces.centroids(:,3));
else
    maxDepth = max(modelRef.G.faces.centroids(:,2));
    minDepth = min(modelRef.G.faces.centroids(:,2));
end
[bc, bd_cells, bd_faces, well_cells, entire_boundary] = introduceBoundaryConditionsFromSchdule(modelRef.G, statesRef, scheduleRef, maxDepth, minDepth, 'wellsPlusSeismic',20);
bd_faces = entire_boundary;
bd_cells = unique(sum(G.faces.neighbors(bd_faces, :), 2));
scheduleRateBhp = flipWellControlsAndAddBC(scheduleRef, statesRef, bc);
% state0.wellSol = statesRef{1}.wellSol;
[wsRateBhp, statesRateBhp] = simulateScheduleAD(state0, modelBhpRate, scheduleRateBhp);

summary_plots2 = plotWellSols({wsRef, wsBhpRate} ,{scheduleRef.step.val, scheduleBhpRate.step.val},...
                            'datasetnames',{'ref. model','rate-bhp model'});
summary_plots3 = plotWellSols({wsBhpRate, wsRateBhp} ,{scheduleBhpRate.step.val, scheduleRateBhp.step.val},...
                            'datasetnames',{'bhp-rate model','rate-bhp model (compatibility)'});

%% Prepare parameters with corresponding bounds 
params_modelBhpRate = [];
params_modelRateBhp = [];

%% Initialize a model with an estimated porosity: here we consider max(pv) from the reference model
InitRock= modelRef.rock;
[maxpv,a] = max(modelRef.operators.pv);

InitRock.poro = InitRock.poro.*0.0 + mean(modelRef.rock.poro);
InitModel = GenericBlackOilModel(modelRef.G, InitRock, modelRef.fluid, 'gas', false);
InitModel.OutputStateFunctions(end+1) = {'ComponentTotalFlux'};
well_data = modelRef.operators.pv(well_cells);
InitModel.operators.pv(well_cells,:)= well_data;
% InitModel.toleranceCNV = 1e-8;
% InitModel.nonlinearTolerance=1.0e-8;
% InitModel.FacilityModel.toleranceWellBHP = 9000;
% InitModel.FacilityModel.toleranceWellRate = 1.0e-8;
% InitModel.FacilityModel.nonlinearTolerance =1.0e-8;
% We run first two simulations with compatible schedule_1 and schedule_2
% with initial parameters
modelBhpRate = InitModel;
modelRateBhp =  InitModel;
[wsBhpRate, statesBhpRate] = simulateScheduleAD(state0, modelBhpRate, scheduleBhpRate);
[wsRateBhp, statesRateBhp] = simulateScheduleAD(state0, modelRateBhp, scheduleRateBhp);

summary_plots3 = plotWellSols({wsBhpRate, wsRateBhp} ,{scheduleBhpRate.step.val, scheduleRateBhp.step.val},...
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
 weighting =  {'WaterRateWeight',  4.4270e+03, ...
'OilRateWeight',    4.4365e+03,...
 'BHPWeight',        (240*barsa)^-1};

%   weighting =  {'WaterRateWeight',  [], ...
%  'OilRateWeight',    [],...
%  'BHPWeight',        []};

initBFGS = true;
itBFGS = 20;
maxIt = 40;
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
                                 'maxIt', itBFGS, 'lbfgsStrategy', 'dynamic', 'lbfgsNum', 20);

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

% compare reference, initial coarse and optimized coarse model outputs
plotWellSols({wsRef, wellSols_opt, wellSols_opt2}, ...
              repmat({scheduleRef.step.val}, 1, 3), ...
              'datasetnames', {'reference','first tuned (Q-N)', 'second tuned (L-M)'});
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