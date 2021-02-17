%% Illustrate use Adjoin optimization for history match
% In this example we explain how to history match two cosared models. One 
% model that have been obtained by upscailing an fine scaled model. This
% upscaled coarse model has production rates very similar to the one
% obtained by the fine scale model however the diference can be minimize by
% a history match. 

% A second coarse model with production result diferent from the fine
% scaled model is also used to show how the robust is the adjoint method
% for history matching

% We history match water and oil production rates and we define as parameters:
%   -Transmisibility of each internal faces, 
%   -Pore volume of each cell
%   -Well index for each well 

mrstModule add agglom upscaling coarsegrid...
               ad-core ad-blackoil ad-props...
               optimization;

use_trans = true;

%% Setting up the fine-scale model
% We make a small model that consists of two different facies with
% contrasting petrophysical properties. 
% Two injectors and two producers are
% placed in the corner at diferent depths.

[model_fine_scale, W, state0] = simpleModelForHMExample();
 model_fine_scale.gas=false;
 model_fine_scale.OutputStateFunctions = {};
 model_fine_scale = model_fine_scale.validateModel();


figure(1), movegui('northwest'); clf,
set(gcf,'Position', [860 450 840 400],'PaperPositionMode','auto');
subplot(1,2,1);
plotCellData(model_fine_scale.G,model_fine_scale.rock.poro,'EdgeAlpha',.5); view(3);
title('Fine-scale model porosity')
plotWell(model_fine_scale.G,W,'Color','k'); axis off tight
pax = [min(model_fine_scale.rock.poro) max(model_fine_scale.rock.poro)];
[hc,hh] = colorbarHist(model_fine_scale.rock.poro, pax,'West');
pos = get(hh,'Position'); set(hh,'Position', pos - [.01 0 0 0]);
pos = get(hc,'Position'); set(hc,'Position', pos - [.03 0 0 0]);
set(gca,'Position',[.12 .075 .315 .9])
drawnow

%% Run simulation for fine-scale model
dt = [1 ,1 ,2 5 10 , 10*ones(1,10)]*day();
schedule = simpleSchedule(dt, 'W', W);

problem = packSimulationProblem(state0, model_fine_scale, schedule, 'model_fine_scale');
[ok, status] = simulatePackedProblem(problem);
[wellSols_fine_scale, states_fine_scale] = getPackedSimulatorOutput(problem);


%% Coarse-scale model
% We make a 3x3x1 coarse grid. 
% TODO: We are upscailing permeability make the same with transmisibility 
% and update this comments:

    % Transmissibilities and well indices are upscaled
    % using two slightly different methods: for the transmissibilities we use
    % global generic boundary conditions and on each coarse face use the
    % solution that has the largest flux orthogonal to the face to compute the
    % upscaled transmissibility. For the wells, we use specific well
    % conditions and use least squares for the flux.
    
    
G = model_fine_scale.G;
rock =  model_fine_scale.rock;
fluid = model_fine_scale.fluid;

hT = computeTrans(G, rock);

p  = partitionCartGrid(G.cartDims, [3 3 1]);
CG = coarsenGeometry(generateCoarseGrid(G, p));

crock = convertRock2Coarse(G, CG, rock);
crock.perm = upscalePerm(G, CG, rock);
[~,~,WC]   = upscaleTrans(CG, hT, 'match_method', 'lsq_flux', ...
                          'bc_method', 'wells_simple', 'wells', W);                     

figure(1), movegui('northwest');
subplot(1,2,2); cla

plotCellData(CG,crock.poro,'EdgeColor','none');
title('Coarse-scale model porosity')

plotFaces(CG, boundaryFaces(CG),...
   'EdgeColor', [0.4 0.4 0.4],'EdgeAlpha',.5, 'FaceColor', 'none'); view(3);
plotWell(G,W,'Color','k'); axis off tight
[hc,hh]=colorbarHist(crock.poro, pax,'East');
pos = get(hh,'Position'); set(hh,'Position', pos + [.01 0 0 0]);
pos = get(hc,'Position'); set(hc,'Position', pos + [.02 0 0 0]);
set(gca,'Position',[.56 .075 .315 .9])
drawnow

 %% Simulating on coarse-scale model
 % After upscailing we create a coarse model and we compare the production
 % results between the fine scaled and coarse scale models
 
%fluid_2 = fluid;
%fluid_2 .krPts  = struct('w', [0 0 1 1], 'ow', [0 0 1 1]);
%scaling = {'SWL', .1, 'SWCR', .2, 'SWU', .9, 'SOWCR', .1, 'KRW', .9, 'KRO', .8};
%scaling = {'SWL', .2, 'SWCR', .2, 'SWU', .8, 'SOWCR', .2, 'KRW', 1, 'KRO', 1};

model_coarse_scale = GenericBlackOilModel(CG, crock, fluid);
model_coarse_scale.gas=false;
%model_coarse_scale = imposeRelpermScaling(model_coarse_scale, scaling{:});

model_coarse_scale.OutputStateFunctions = {};
model_coarse_scale = model_coarse_scale.validateModel();
 
s0 = [0,1];
p0 = 100*barsa;
state0  = initState(CG, WC,p0,s0);

schedule = simpleSchedule(dt, 'W', WC);
[wellSols_coarse_scale,states_coarse_scale] = simulateScheduleAD(state0, model_coarse_scale, schedule);

summary_plots = plotWellSols({wellSols_fine_scale,wellSols_coarse_scale},{schedule.step.val,schedule.step.val});
movegui('northeast')
drawnow
%% Preparing parameters and scaling values for each one
% We define each parameters well index transmisibility, pore volume, and
% it's scailing values.

prob = struct('model', model_coarse_scale, 'schedule', schedule, 'state0', state0);

parameters{1} = ModelParameter(prob, 'name', 'conntrans','relativeLimits', [.01 1.5]);
parameters{2} = ModelParameter(prob, 'name', 'porevolume','relativeLimits', [.01 3]);
parameters{3} = ModelParameter(prob, 'name', 'transmissibility','relativeLimits', [.1 2]);
%% 
values = applyFunction(@(p)p.getParameterValue(prob), parameters);
% scale values
u = cell(size(values));
for k = 1:numel(u)
    u{k} = parameters{k}.scale(values{k});
end
p0_ups = vertcat(u{:});  
  
% Defining the weights to evaluate the match
weighting =  {'WaterRateWeight',  (20/day)^-1, ...
              'OilRateWeight',    (20/day)^-1, ...
              'BHPWeight',        (100*barsa)^-1};            

 obj = @(model, states, schedule, states_ref, tt, tstep, state) matchObservedOW(model, states, schedule, states_fine_scale,...
           'computePartials', tt, 'tstep', tstep, weighting{:},'state',state,'from_states',false);

 objScaling = 1;       
 [misfitVal_0,gradient,wellSols_0,states_0] = evaluateMatch_simple(p0_ups,obj,state0,model_coarse_scale,schedule,objScaling,parameters, states_fine_scale);

                                                           

  
obj_scaling     = abs(misfitVal_0);      % objective scaling  
objh = @(p)evaluateMatch_simple(p, obj, state0, model_coarse_scale, schedule, obj_scaling ,parameters,  states_fine_scale);

figure(10).reset; movegui('south');
[v, p_opt, history] = unitBoxBFGS(p0_ups, objh,'gradTol',             1e-2, ...
                                              'objChangeTol',        0.5e-3,...
                                              'maxIt',               8);
[misfitVal_opt,gradient_opt,wellSols_opt] = evaluateMatch_simple(p_opt, obj, state0, model_coarse_scale, schedule, obj_scaling ,parameters, states_fine_scale);
 
figure(summary_plots.Number)
plotWellSols({wellSols_fine_scale,wellSols_0,wellSols_opt},...
              {schedule.step.val,schedule.step.val,schedule.step.val},...
              'datasetnames',{'fine scale model','initial upscaled model','history matched upscaled model'},...
              'linestyles', {'o', '--', '-'},...
              'figure',summary_plots.Number)
drawnow
 %% Optimization 2:  initial upscaled model with a ramdom perturbation
 
 
% Adding a perturmation to the initial values
  for k = 1:numel(u)
    values{k} =  values{k}+1.5*(rand(size(values{k}))-0.5).*values{k};
    prob = parameters{k}.setParameterValue(prob,values{k});
    u{k} = parameters{k}.scale(values{k});
  end
  p0_mean = vertcat(u{:});
  model_coarse_scale = prob.model;
  schedule = prob.schedule;
  state0   = prob.state0;
 
       
 [misfitVal_0,gradient,wellSols_0_mean,states_0] = evaluateMatch_simple(p0_mean,obj,state0,model_coarse_scale,schedule,objScaling,parameters, states_fine_scale);

obj_scaling     = abs(misfitVal_0);      % objective scaling  
objh = @(p)evaluateMatch_simple(p,obj, state0, model_coarse_scale, schedule, obj_scaling ,parameters,  states_fine_scale);

gco(figure(10))
figure(10).reset, movegui('south')
[v, p_opt, history] = unitBoxBFGS(p0_mean, objh,'gradTol',            1e-2, ...
                                              'objChangeTol',        0.5e-3,...
                                              'maxIt',               4);
[misfitVal_opt,gradient_opt,wellSols_opt_mean] = evaluateMatch_simple(p_opt, obj, state0, model_coarse_scale, schedule, obj_scaling ,parameters, states_fine_scale);

%clf(summary_plots)
figure(summary_plots.Number)
% gcbf(figure(summary_plots.Number))
% summary_plots
plotWellSols({wellSols_fine_scale,wellSols_0,wellSols_opt,wellSols_0_mean,wellSols_opt_mean},...
              {schedule.step.val,schedule.step.val,schedule.step.val,schedule.step.val,schedule.step.val},...
              'datasetnames',{'fine scale model','initial upscaled model','history matched upscaled model','perturbed upscaled model','history matched perturbed upscaled model'},...
              'linestyles', {'o', '--', '-','--', '-'},...
              'figure',summary_plots.Number)
