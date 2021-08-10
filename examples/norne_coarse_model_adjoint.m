mrstModule add ad-core ad-blackoil deckformat ...
               agglom upscaling coarsegrid...
               mrst-gui ad-props incomp optimization...
               network-models example-suite linearsolvers 

%% Setting up the reference model           
setting_up_simplified_norne_model

wellSols_ref =  wellSols;
model_ref    = example.model;
states_ref   = states;
state0_ref   = example.state0;
schedule_ref = example.schedule;
W_ref        = schedule_ref.control.W;


figure(1), clf, subplot(1,2,1)
plotCellData(model_ref.G,model_ref.rock.poro,'EdgeAlpha',.5); view(3);
title('Fine-scale model porosity')
plotWell(model_ref.G,W_ref,'Color','k'); axis off tight
drawnow

time_steps     = schedule_ref.step.val(1:6); 
training_steps = 1:6;
partition      = [6 8 1];

%% Coarse-scale model
% We make a coarse grid defined by partition, and perform a simple
% upscaling to obtain a coarse model
%p  = partitionCartGrid(model_ref.G.cartDims, partition );
p = partitionUI(model_ref.G,partition);
model_c = upscaleModelTPFA(model_ref, p);
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
schedule_training   = simpleSchedule(time_steps(training_steps), 'W', W_ref);
schedule_training_c = upscaleSchedule(model_c, schedule_training);


figure(1)
subplot(1,2,2); cla;
plotCellData(model_c.G, model_c.rock.poro,'EdgeColor','none');
title('Coarse-scale model porosity')
plotFaces(model_c.G, boundaryFaces(model_c.G), 'EdgeColor', [0.4 0.4 0.4], ...
         'EdgeAlpha',.5, 'FaceColor', 'none'); view(3);
plotWell(model_ref.G, W_ref, 'Color', 'k'); axis off tight


 %% Simulate initial upscaled coarse model for full time
state0_c   = upscaleState(model_c, model_ref, state0_ref);
schedule_c = upscaleSchedule(model_c, schedule_ref);
[wellSols_c, states_c] = simulateScheduleAD(state0_c, model_c, schedule_c);
summary_plots = plotWellSols({wellSols_ref, wellSols_c} ,{schedule_ref.step.val, schedule_c.step.val},...
                            'datasetnames',{'fine scale model','initial upscaled model'});
drawnow

 %% Preparing parameters
nc =  model_c.G.cells.num;

% select configuration for sensitivity computations
prob_training = struct('model', model_c, 'schedule', schedule_training_c, 'state0', state0_c);
config = {...%name      include     scaling    boxlims     lumping     subset       relativeLimits
          'porevolume',       1,   'linear',       [],  [],   [],   [.001 4]
          'conntrans',        1,   'log',          [],          [],       [] ,      [.001 100]
          'transmissibility', 1,   'log'   ,       [],  [], [] ,   [.001 100]
          'swl',              1,   'linear',       [],  ones(nc,1),       [],       [.5 2]
          'swcr',             1,   'linear',       [],  ones(nc,1),       [],       [.5 2]
          'swu',              1,   'linear',       [],  ones(nc,1),       [],       [.5 2]
          'sowcr',            1,   'linear',       [],  ones(nc,1),       [],       [.5 2]
          'krw',              1,   'linear',       [],  ones(nc,1),       [],       [.5 2]
          'kro',              1,   'linear',       [],  ones(nc,1),       [],       [.5 2]
          'sw',               1,   'linear',   [0 .6],  [],   [],  [0 10]
          'pressure'          1,   'linear',       [],  [],   [],  [.1 4]};
parameters = [];
for k = 1:size(config,1)
    if config{k, 2} ~= 0
        parameters = addParameter(parameters, prob_training, 'name',    config{k,1}, ...
                                                     'scaling', config{k,3}, ...
                                                     'boxLims', config{k,4}, ...
                                                     'lumping', config{k,5}, ...
                                                     'subset',  config{k,6},...
                                                     'relativeLimits',config{k,7});
    end
end
 



 %% Simulation of the reference model and Data-driven model

         weighting =  {'WaterRateWeight',  (500/day)^-1, ...
                       'OilRateWeight',    (500/day)^-1, ...
                       'BHPWeight',        (50*barsa())^-1};
 
p0 = getScaledParameterVector(prob_training, parameters);
 
  obj = @(model, states, schedule, states_ref, tt, tstep, state) matchObservedOW(model, states, schedule, states_ref,...
            'computePartials', tt, 'tstep', tstep, weighting{:},'state',state,'from_states',false);
       
 [misfitVal_0,~,wellSols_0,states_0] = evaluateMatch(p0,obj,prob_training,parameters, states_ref,'Gradient','none');          
 
 plotWellSols({wellSols_ref,wellSols_0},{schedule_ref.step.val,prob_training.schedule.step.val})
  



%% Optimization
  
objh = @(p)evaluateMatch(p,obj,prob_training,parameters,states_ref);

[v, p_opt, history] = unitBoxBFGS(p0, objh,'objChangeTol',  1e-8, 'maxIt',30, 'lbfgsStrategy', 'dynamic', 'lbfgsNum', 5);

%% Simulating all simulation time
prob = prob_training;
prob.schedule = simpleSchedule(schedule_ref.step.val, 'W', prob.schedule.control.W);
 
 [~,~,wellSols_opt] = evaluateMatch(p_opt,obj,prob,parameters, states_ref,'Gradient','none');
 [~,~,wellSols_0] = evaluateMatch(p0,obj,prob,parameters, states_ref,'Gradient','none');


fh = plotWellSols({wellSols_ref,wellSols_0,wellSols_opt},{schedule_ref.step.val,prob.schedule.step.val,prob.schedule.step.val})
set(fh, 'name','Norne')
legend('reference model','initial DD model','optimize DD model')
                                       