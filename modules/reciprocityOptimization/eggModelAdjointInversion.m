clear all;
%gravity reset off;
%% Parameter tuning of a very coarse upscaling of the Egg model   
mrstModule add ad-core ad-blackoil deckformat agglom upscaling coarsegrid...
        mrst-gui ad-props incomp optimization test-suite linearsolvers 
%% Setup reference model
gravity reset off;
test    = TestCase('egg_wo', 'realization', 1);
problem = test.getPackedSimulationProblem();
problem.SimulatorSetup.model.gravity = [0 0 0];
% We focus on the first 60 steps
schedule = problem.SimulatorSetup.schedule;
schedule.step.val     = schedule.step.val(1:60);
schedule.step.control = schedule.step.control(1:60);
Wref     = schedule.control.W;
dt       = schedule.step.val;
nstep    = numel(schedule.step.val);

% For the training run, we create a schedule with varying controls
perturbStep = [ones(1,4), round(.5+(2:(nstep-3))/2)]';
rng(0)
schedule = perturbedSimpleSchedule(dt, 'W', Wref, ...
   'pressureFac', .01, 'rateFac', .4, 'perturbStep', perturbStep);
problem.SimulatorSetup.schedule = schedule;

% Simulate
rockRef = problem.SimulatorSetup.model.rock;
rockRef.poro = rockRef.perm(:,1)./max(rockRef.perm(:,1)).*0.3+0.1;
problem.SimulatorSetup.model.rock= rockRef;
modelRef = problem.SimulatorSetup.model.setupOperators(problem.SimulatorSetup.model.G, rockRef);
problem.SimulatorSetup.model = modelRef;
%simulatePackedProblem(problem,'restartStep',60);
gravity off;
%[wsRef, statesRef] = getPackedSimulatorOutput(problem);
modelRef    = problem.SimulatorSetup.model;

%% Coarse-scale model
% We make a coarse grid defined by a uniform 11 x 11 x 7 partition 
blockIx = partitionUI(modelRef.G, [10, 10, 7]);
blockIx = processPartition(modelRef.G, blockIx);
blockIx = compressPartition(blockIx);
% Perform a simple upscaling to obtain a coarse model
modelCoarse = upscaleModelTPFA(modelRef, blockIx);
%modelCoarse.rock.poro = modelCoarse.rock.perm(:,3)./max(modelCoarse.rock.perm(:,3));
%modelCoarse.rock.poro(modelCoarse.rock.poro>0.75) = 0.75;
%modelCoarse.rock.poro = modelCoarse.rock.perm(:,1)./max(modelCoarse.rock.perm(:,1)).*0.3+0.1;;
%modelCoarse.AutoDiffBackend = AutoDiffBackend();
modelCoarse.OutputStateFunctions{end+1} ='Mobility';
% We want to include rel-perm scaling as tunabale parameters, so include
% these for the coarse model. These parameters have no effect for the
% initial coarse model (they are set equal to the ones given by the
% rel-perm curves).
% pts = modelCoarse.fluid.krPts;
% scaling = {'SWL',   pts.w(1), 'SWCR', pts.w(2), 'SWU', pts.w(3), ...
%            'SOWCR', pts.ow(2), 'KRW',  pts.w(4), 'KRO', pts.ow(4)};
% modelCoarse = imposeRelpermScaling(modelCoarse, scaling{:});
%modelCoarse.toleranceCNV = 1e-10;  % tighter tolerance to improve gradient accuracy
modelCoarse.nonlinearTolerance=1.0e-9;
modelCoarse.toleranceCNV = 1.0e-10;
modelCoarse.toleranceMB = 1.0e-7;
modelCoarse.FacilityModel.toleranceWellBHP = 1;
modelCoarse.FacilityModel.toleranceWellRate = 1.0e-11;
modelCoarse.FacilityModel.nonlinearTolerance =1.0e-9;
%% Plot reference and coarse model grids with wells
figure, subplot(1,2,1)
plotGrid(modelRef.G, 'EdgeAlpha',.2); 
title('Fine-scale grid (18553 cells)')
plotWell(modelRef.G, Wref, 'Color', 'k', 'FontSize', 10);
axis off tight, view(174,60), camlight headlight
subplot(1,2,2)
plotGrid(modelCoarse.G, 'EdgeAlpha',.8);
plotWell(modelRef.G, Wref, 'Color', 'k', 'FontSize', 10);
axis off tight, view(174,60), camlight headlight
title('Coarse-scale grid (33 cells)')
%% Simulate initial upscaled coarse model and compare to reference
%stateCoarse0   = upscaleState(modelCoarse, modelRef, test.state0);
stateCoarse0   = upscaleState(modelCoarse, modelRef, test.state0);
%stateCoarse0.pressure   = zeros(modelCoarse.G.cells.num,1)  + max(test.state0.pressure);
stateCoarse0.s   = zeros(modelCoarse.G.cells.num,2)  + test.state0.s(1,:);
scheduleCoarse = upscaleSchedule(modelCoarse, schedule,'wellUpscaleMethod', 'sum');%% Specify parameters for tuning
%modelCoarse.rock.poro = modelCoarse.rock.perm(:,1)./max(modelCoarse.rock.perm(:,1)).*0.3+0.1;
%modelCoarse = modelCoarse.setupOperators(modelCoarse.G, modelCoarse.rock);
%modelCoarse = GenericBlackOilModel(modelCoarse.G, modelCoarse.rock, modelCoarse.fluid, 'gas', false);
setupCoarse = struct('model', modelCoarse, 'schedule', scheduleCoarse, 'state0', stateCoarse0);
pv = setupCoarse.model.operators.pv;
% set up 'matrix' for parameter options for easier editing. The specific
% limits set for the various parameters influences the tuning/optimization 
% procedure to a large extent
[wsCoarse, statesCoarse] = simulateScheduleAD(setupCoarse.state0, setupCoarse.model, setupCoarse.schedule);
stateCoarse0 = setupCoarse.state0;
scheduleCoarse = setupCoarse.schedule;
modelCoarse  = setupCoarse.model;
fluidCoarse  = setupCoarse.model.fluid;
%plotWellSols({wsRef, wsCoarse}, ...
%    {schedule.step.val, scheduleCoarse.step.val},...
%    'datasetnames',{'fine scale model','initial upscaled model'});
%gravity reset off;
 %% Simulate initial upscaled coarse model for full time
state0_1   = stateCoarse0;
schedule_1= scheduleCoarse;
model_1 = setupCoarse.model;
[wellSols_1, states_1] = simulateScheduleAD(state0_1,model_1 , schedule_1);

% create a second schedule
WW=scheduleCoarse.control.W;
for i =1:length(WW)
    if strcmp(WW(i).type,'rate')
        WW(i).type     = 'bhp';    
        WW(i).sign     = 1;
        WW(i).lims  = []; 
    else       
        WW(i).type     = 'rate';
        WW(i).lims  = [];
        WW(i).sign     = -1;
    end
end
cells=[];
for j = 1:length(schedule_1.control(1).W)
    cells=[cells;schedule_1.control(1).W(j).cells];
end
bc =[];
G = modelCoarse.G;
[neighborship, n_isnnc] = getNeighbourship(G, 'Topological', true);
[cellNo, cf, cn_isnnc] = getCellNoFaces(G);
nif    = size(neighborship, 1);
ncf    = size(cf, 1);
nc     = G.cells.num;
nw     = length(WW);
n      = nc + nw;
bdf=[find(abs(G.faces.centroids(:,3) - 4000)<1.01);find(abs(G.faces.centroids(:,3) -4028)<1.21)];
cells_bdf= unique(sum(G.faces.neighbors(bdf,:),2));
[a,b]=intersect(cells_bdf,cells);
time_steps = scheduleCoarse.step.val;
[rest, rest_i]=  setdiff(cells_bdf,a);
random_indices =randperm(length(rest), 24);
%b = 1:length(cells);
b = [b;rest_i(random_indices)];
for i =1:length(time_steps)
     bc = addBC([], bdf(b), 'pressure',        ...
            statesCoarse{i}.pressure(cells_bdf(b)),'sat', [statesCoarse{i}.s(cells_bdf(b),:)]);    
    
WW=scheduleCoarse.control(scheduleCoarse.step.control(i)).W;
for ii =1:length(WW)
    if strcmp(WW(ii).type,'rate')
        WW(ii).type     = 'bhp';    
       % WW(ii).sign     = 1;
        WW(ii).lims  = []; 
    else       
        WW(ii).type     = 'rate';
        WW(ii).lims  = [];
       % WW(ii).sign     = -1;
    end
end

    for j =1:length(WW)    
        if strcmp(WW(j).type,'rate')
            WW(j).val      = statesCoarse{i}.wellSol(j).qOs + statesCoarse{i}.wellSol(j).qWs;            
        else           
            WW(j).val = statesCoarse{i}.wellSol(j).bhp;
        end

    end
    schedule2{i} = simpleSchedule(time_steps(i), 'W', WW,'bc',bc); 
end
schedule_2 = combineSchedules(schedule2{:}, 'makeConsistent', false);
[wellSols_2, states_2] = simulateScheduleAD(state0_1, model_1, schedule_2);


%% Prepare parameters with corresponding bounds 
params_model1 = [];
params_model2 = [];
poro = modelCoarse.rock.poro;
% Initial permeability
rock_f_I = modelCoarse.rock;
rock_f_I.poro = modelCoarse.rock.poro.*0 +  max(poro);
rock_f_I.poro(cells) = modelCoarse.rock.poro(cells);
model_f_I = GenericBlackOilModel(modelCoarse.G, rock_f_I, modelCoarse.fluid, 'gas', false);
% model_f_I.toleranceCNV = 1e-7;
% model_f_I.nonlinearTolerance=1.0e-7;
% %model_f_I.FacilityModel.toleranceWellBHP = 90000;
% model_f_I.FacilityModel.toleranceWellRate = 1.0e-7;
% model_f_I.FacilityModel.nonlinearTolerance =1.0e-7;
% We run first two simulations with compatible schedule_1 and schedule_2 
model_1 = model_f_I;
model_2 = model_f_I;
[wellSols_1, states_1] = simulateScheduleAD(state0_1, model_1, schedule_1);
[wellSols_2, states_2] = simulateScheduleAD(state0_1, model_2, schedule_2);

setup_model1 = struct('model', model_f_I, 'schedule', schedule_1, 'state0', state0_1);
setup_model2 = struct('model', model_f_I, 'schedule', schedule_2, 'state0', state0_1);
% a = min(modelCoarse.rock.poro(:,1))./(rock_f_I.poro);
% b = max(modelCoarse.rock.poro(:,1))./(rock_f_I.poro);
a = min(modelCoarse.operators.pv)-0.01;
b = max(model_f_I.operators.pv)+0.01;
params_model1 = addParameter(params_model1,setup_model1, 'name', 'porevolume','boxLims', [a b]);
% params_model1 = addParameter(params_model1,setup_model1, 'name', 'transmissibility', 'scaling','log','boxLims', [0.18 0.82]);
params_model2 = addParameter(params_model2,setup_model2, 'name', 'porevolume','boxLims', [a b]);
% params_model2 = addParameter(params_model2,setup_model2, 'name', 'transmissibility', 'scaling', 'log','boxLims', [0.18 0.82]);

%% Coarse model tuning
p0 = getScaledParameterVector(setup_model1, params_model1);
% % Weight production rates/bhp relative to their magnitudes
% weighting =  {'WaterRateWeight',  (4*stb/day)^-1, ...
%               'OilRateWeight',    (4*stb/day)^-1,...
%               'BHPWeight',        (45*barsa)^-1};   
weighting  = {'WaterRateWeight',  [], ...bc
              'OilRateWeight',    [], ...
              'BHPWeight',        []};
% define objective
% method =2;
% if (method ==1)
obj = @(model_1, states_1, schedule_1, model_2, states_2, schedule_2,states_c, doPartials, tstep, state1,state2)...
    matchObservedTwinModelsOW_saved(model_1, states_1, schedule_1, model_2, states_2, schedule_2,states_c,cells_bdf,bdf,cells,modelCoarse.operators.pv(cells,:), ...
    'computePartials', doPartials, 'tstep', tstep, 'state1', state1,'state2', state2, ...
    'from_states', false, weighting{:});
% provide handle to function that runs simulation for provided parameters, 
% computes objective and gradient (if requested).
objh = @(p)evaluateMatchTwinModels(p, obj, setup_model1,setup_model2, params_model1, statesCoarse,cells,modelCoarse.operators.pv(cells,:));
% run optimization
[v, p_opt, history] = unitBoxBFGS(p0, objh, 'objChangeTol', 1e-7, 'gradTol', 1e-7, ...
                                 'maxIt', 10, 'lbfgsStrategy', 'dynamic', 'lbfgsNum', 20);

%% Re-run simulation for optimal parameters for full time-horizon
setup_opt = updateSetupFromScaledParameters(setup_model1, params_model1, p_opt);
setup_opt.schedule.step = scheduleCoarse.step; % full time-horizon
[wellSols_opt, states_opt] = simulateScheduleAD(setup_opt.state0, setup_opt.model, setup_opt.schedule);
% else (method ==2)
%%%%%%%%%%%
%% Model calibration Levenberg-Marquard (using full Jacobian)
mismatchFn2 = @(model_1, states_1, schedule_1, model_2, states_2, schedule_2,states_c, doPartials, tstep, state1,state2)...
    matchObservedTwinModelsOW(model_1, states_1, schedule_1, model_2, states_2, schedule_2,states_c,cells_bdf,bdf,cells,modelCoarse.operators.pv(cells,:), ...
    'computePartials', doPartials, 'tstep', tstep, 'state1', state1,'state2', state2, ...
    'from_states', false, weighting{:}, 'mismatchSum', false);
objh2 = @(p) evaluateMatchTwinModelsSummands(p, mismatchFn2, setup_model1,setup_model2, params_model1, statesCoarse,cells,modelCoarse.operators.pv(cells,:));
% The calibration can be improved by taking a large number of iterations,
% but here we set a limit of 30 iterations
[v2, p_opt2, history2] = unitBoxLM(p_opt, objh2, 'maxIt', 100);
% end



% %% Specify parameters for tuning
% setup_init = struct('model', modelCoarse, 'schedule', scheduleCoarse, 'state0', stateCoarse0);
% pv = modelCoarse.operators.pv;
% % set up 'matrix' for parameter options for easier editing. The specific
% % limits set for the various parameters influences the tuning/optimization 
% % procedure to a large extent
% config = {...
%      %name           include    scaling              boxlims   relativeLimits  
%     'porevolume',       1,     'linear',    [.01*pv, 1.5*pv],              []   
%     'conntrans',        1,        'log',                  [],     [1e-2, 1e2]      
%     'transmissibility', 1,        'log',                  [],     [1e-2  1e2]  
%     'swl',              1,     'linear',             [0, .3],              []
%     'swcr',             1,     'linear',             [0, .4],              []
%     'swu',              1,     'linear',             [.7, 1],              []
%     'sowcr',            1,     'linear',             [0, .4],              []
%     'krw',              1,     'linear',           [.5, 1.5],              []
%     'kro',              1,     'linear',           [.5, 1.5],              []};
% parameters = [];
% for k = 1:size(config,1)
%     if config{k, 2} == 0, continue, end
%     parameters = addParameter(parameters, setup_init, ...
%         'name',    config{k,1}, 'scaling', config{k,3}, ...
%         'boxLims', config{k,4}, 'relativeLimits',config{k,5});
% end
% 
% 
% 
% %% Define the mismatch function
% % Function weighting influences the match of each quantity. Rate-weighting
% % should be on the same order as (inverse of) rates. BHP-weighting on the
% % order of pressure drop in the model.
% weighting  = {'WaterRateWeight',  1/(150/day), ...
%               'OilRateWeight',    1/(80/day), ...
%               'BHPWeight',        1/(20*barsa)};
% % make handle          
% mismatchFn = @(model, states, schedule, states_ref, compDer, tstep, state) ...
%     matchObservedOW(model, states, schedule, states_ref,...
%                    'computePartials', compDer, 'tstep', tstep, weighting{:},...
%                    'state', state, 'from_states', false);
% 
% %% Model calibration Quasi-Newton
% pvec = getScaledParameterVector(setup_init, parameters);
% objh = @(p) evaluateMatch(p, mismatchFn, setup_init, parameters, statesRef);
% % The calibration can be improved by taking a large number of iterations,
% % but here we set a limit of 30 iterations
% [v1, p_opt1, history1] = unitBoxBFGS(pvec, objh, 'objChangeTol', 1e-5, ...
%                                   'maxIt', 20, 'logPlot', true);
% 
% %% Model calibration Levenberg-Marquard (using full Jacobian)
% mismatchFn2 = @(model, states, schedule, states_ref, compDer, tstep, state) ...
%     matchObservedOW(model, states, schedule, states_ref, ...
%         'computePartials', compDer, 'tstep', tstep, weighting{:},...
%         'state', state, 'from_states', false, 'mismatchSum', false);
% objh2 = @(p) evaluateMatchSummands(p, mismatchFn2, setup_init, parameters, statesRef);
% % The calibration can be improved by taking a large number of iterations,
% % but here we set a limit of 30 iterations
% [v2, p_opt2, history2] = unitBoxLM(pvec, objh2, 'maxIt', 20);
% 
% %% Compare convergence history
% figure
% semilogy(abs(history1.val(1:20))', '-o', 'LineWidth', 2);
% hold on
% semilogy(abs(history2.val(1:20))', '-o', 'LineWidth', 2);
% set(gca, 'Fontsize', 12); grid on
% xlabel('Iteration'), ylabel('Residual')
% legend({'Quasi-Newton', 'Levenberg-Marquard'})
% 
% %% Create new coarse model setups with the optimized parameters, and rerun 
% %  for the optimized parameters
% setup_opt1 = updateSetupFromScaledParameters(setup_init, parameters, p_opt1); 
% [wellSols_opt1, states_opt1] = simulateScheduleAD(setup_opt1.state0, setup_opt1.model, setup_opt1.schedule);
% setup_opt2 = updateSetupFromScaledParameters(setup_init, parameters, p_opt2); 
% [wellSols_opt2, states_opt2] = simulateScheduleAD(setup_opt2.state0, setup_opt2.model, setup_opt2.schedule);
% % compare reference, initial coarse and optimized coarse model outputs
% plotWellSols({wsRef, wsCoarse, wellSols_opt1, wellSols_opt2}, ...
%               repmat({schedule.step.val}, 1, 4), ...
%               'datasetnames', {'reference','coarse initial','coarse tuned (Q-N)', 'coarse tuned (L-M)'});
% 
% %% Plot the pore volume updates
% % fetch pore volume differences in initial and tuned coarse models
% dpv = setup_opt2.model.operators.pv - setup_init.model.operators.pv;
% figure
% plotCellData(modelCoarse.G, dpv, 'EdgeColor','none');
% plotFaces(modelCoarse.G, boundaryFaces(modelCoarse.G), 'EdgeColor', [0.4 0.4 0.4], ...
%          'EdgeAlpha',.5, 'FaceColor', 'none');
% view(174,60);
% plotWell(modelRef.G, Wref, 'Color', 'k', 'FontSize', 10); axis off tight
% colorbar('south');
%                         
% %% Compare reference, initial coarse and optimizes coarse for a different schedule 
% rng(100);
% W = test.schedule.control.W;
% s_new = perturbedSimpleSchedule(dt, 'W', Wref, ...
%     'pressureFac', .01, 'rateFac', .2, 'perturbStep', ones(numel(dt),1));
% ws_new1 = simulateScheduleAD(problem.SimulatorSetup.state0, modelRef, s_new, ...
%     'NonLinearSolver', problem.SimulatorSetup.NonLinearSolver);
% 
% w_opt = setup_opt2.schedule.control(1).W;
% w_old = scheduleCoarse.control(1).W;   
% s_opt = simpleSchedule(dt, 'W', w_opt);
% s_old = simpleSchedule(dt, 'W', w_old);
% for kw = 1:numel(Wref)
%     s_opt.control.W(kw).val = s_new.control.W(kw).val;
%     s_old.control.W(kw).val = s_new.control.W(kw).val;
% end
% ws_coarse_new1 = simulateScheduleAD(setup_opt2.state0, setup_opt2.model, s_opt);
% ws_coarse_old1 = simulateScheduleAD(stateCoarse0, modelCoarse, s_old);
% plotWellSols({ws_new1,  ws_coarse_old1, ws_coarse_new1}, ...
%     repmat({schedule.step.val}, 1, 3), ...
%     'datasetnames', {'reference','coarse initial','coarse tuned'});
% 
% %% Copyright Notice
% %
% % <html>
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
