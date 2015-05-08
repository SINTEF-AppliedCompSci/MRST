function [optim, init, history] = optimizeRates(initState, model, schedule, ...
                                                min_rates, max_rates, varargin)
   
% optim.obj_val_steps
% optim.obj_val_total
% optim.schedule
% optim.wellSols
% optim.states
   
   moduleCheck('ad-core', 'ad-props', 'optimization');

   opt.obj_scaling = [];  % Compute explicitly if not provided from outside
   opt.leak_penalty = 10; % Only used if 'obj_fun' not provided 
   opt.obj_fun = @(wellSols, states, schedule, varargin) ...
                  leak_penalizer(model, wellSols, states, schedule, opt.leak_penalty, ...
                                 varargin{:});
   
   opt = merge_options(opt, varargin{:});
   
   num_wells = numel(schedule.control(1).W);
   assert(numel(max_rates) == num_wells);
   assert(numel(min_rates) == num_wells);
   
   init.schedule = schedule;

   %% Compute initial objective value (if required for scaling)
   if isempty(opt.obj_scaling)
      [init.wellSols, init.states] = ...
          simulateScheduleAD(initState, model, schedule); 
      
      init.obj_val_steps = cell2mat(opt.obj_fun(init.wellSols, ...
                                                init.states, ...
                                                init.schedule));
      init.obj_val_total = sum(init.obj_val_steps);
      
      % Use value of objective function before optimization as scaling.
      opt.obj_scaling = init.obj_val_total;
   end
   
   
   %% Define limits, scaling and objective function
   
   scaling.boxLims = [min_rates(:), max_rates(:)];
   scaling.obj     = opt.obj_scaling;
   
   obj_evaluator = @(u) evaluate_objective(u, opt.obj_fun, model, initState, ...
                                           schedule, scaling); 
   
   %% Call optimization routine
   
   u = schedule2control(schedule, scaling);
   u = u(1:end-num_wells); % ignore last well configuration, which represent
                           % migration stage (wells should be shut off)
   [opt_val, u_opt, history] = unitBoxBFGS(u, obj_evaluator);
   
   %% Preparing solution structures
   
   optim.schedule = control2schedule([u_opt; zeros(num_wells,1)], schedule, scaling);

   [optim.wellSols, optim.states] = simulateScheduleAD(initState, ...
                                                       model, ...
                                                       optim.schedule);
   
   optim.obj_val_steps = opt.obj_fun(optim.wellSols, optim.states, optim.schedule);
   optim.obj_val_total = sum(cell2mat(optim.obj_val_steps));

   assert(opt_val == optim.obj_val_total);
end

% ----------------------------------------------------------------------------

function obj = leak_penalizer(model, wellSols, states, schedule, leak_penalty, varargin)
% obj is here a cell array with one cell per timestep in the schedule
   
   opt.ComputePartials = false;
   opt.tStep = [];
   opt = merge_options(opt, varargin{:});

   assert(xor(opt.ComputePartials, isempty(opt.tStep)));
   
   pvol = model.G.cells.volumes .* model.G.cells.H .* model.rock.poro;
   if isfield(model.fluid, 'pvMultR')
      pvMult = model.fluid.pvMultR;
   else
      pvMult = @(p) ones(size(p));
   end
   
   numSteps = numel(schedule.step.val);
   assert(numel(states) == numSteps);
   obj = repmat({[]}, numSteps, 1);

   total_injected = zeros(numSteps, 1);
   total_leaked = zeros(numSteps, 1);
   prev_injected = 0;
   prev_total_leaked = 0;
   for step = 1:numSteps
      st = states{step};
      dt = schedule.step.val(step);
      
      inx = ([wellSols{step}.sign] > 0)';
      qGs = [wellSols{step}.qGs]';
      P   = st.pressure;
      sG  = st.s(:,2);

      recently_injected = sum(qGs(inx)) * dt;
      total_injected(step) = prev_injected + recently_injected;
      prev_injected = total_injected(step);
      
      
      total_in_place = sum(pvol .* pvMult(P) .* sG .* model.fluid.bG(P));

      total_leaked(step) = total_injected(step) - total_in_place;
      
      recently_leaked = total_leaked(step) - prev_total_leaked;
      prev_total_leaked = total_leaked(step);
                        
      % Objective value for this timestep equals the amount injected during
      % this step, less the amount leaked during the step multiplied by a
      % weight. 
      %obj{step} = recently_injected - leak_penalty * recently_leaked;
      obj{step} = -recently_injected;
   end
   
   if ~isempty(opt.tStep)
      assert(opt.ComputePartials);

      st = states{opt.tStep};
      dt = schedule.step.val(opt.tStep);
      
      inx = ([wellSols{opt.tStep}.sign] > 0)';
      qGs = [wellSols{opt.tStep}.qGs]';
      P   = st.pressure;
      sG  = st.s(:,2);
      
      if opt.ComputePartials
         numWells = numel(qGs);
         bhp = zeros(numWells, 1); % placeholder
         qWs = zeros(numWells, 1); % placeholder
         [P, sG, dummy, dummy, qGs] = initVariablesADI(P, sG, bhp, qWs, qGs); %#ok
      end

      sum_qGs = qGs(inx);
      if numel(double(sum_qGs)) > 1
         % Hack to get around a bug in ADI @@ report this!
         sum_qGs = sum(sum_qGs);
      end
      
      recently_injected = sum_qGs * dt;
      
      total_in_place = sum(pvol .* pvMult(P) .* sG .* model.fluid.bG(P));
      
      tot_inj = recently_injected;
      if opt.tStep > 1
         tot_inj = tot_inj + total_injected(opt.tStep-1);
      end
      
      recently_leaked = tot_inj - total_in_place;
      
      if opt.tStep > 1
         recently_leaked = recently_leaked - total_leaked(opt.tStep-1);
      end
      
      %obj = recently_injected - leak_penalty * recently_leaked;
      obj = -recently_injected;
   end
   
end


% ----------------------------------------------------------------------------

function [val, der, wellSols, states] = ...
       evaluate_objective(u, obj_fun, model, initState, schedule, scaling) 
   
   minu = min(u);
   maxu = max(u);
   if or(minu < -eps , maxu > 1+eps)
      warning('Controls are expected to lie in [0 1]\n')
   end

   boxLims = scaling.boxLims;
   if isfield(scaling, 'obj')
      objScaling = scaling.obj;
   else
      objScaling = 1;
   end
   
   % update schedule:
   Wnum = numel(schedule.control(1).W);
   schedule = control2schedule([u; zeros(Wnum,1)], schedule, scaling);
   
   % run simulation:
   [wellSols, states] = simulateScheduleAD(initState, model, schedule);
   
   % compute objective:
   vals = obj_fun(wellSols, states, schedule);
   val  = sum(cell2mat(vals))/abs(objScaling);

   % run adjoint:
   if nargout > 1
      objh = @(tstep)obj_fun(wellSols, states, schedule, 'ComputePartials', true, 'tStep', tstep);
      g    = computeGradientAdjointAD(initState, states, model, schedule, objh);
      % scale gradient:
      der = scaleGradient(g, schedule, boxLims, objScaling);
      der = der(1:end - Wnum);
      der = vertcat(der{:});
   end
end

% ----------------------------------------------------------------------------

function grd = scaleGradient(grd, schedule, boxLims, objScaling)
   dBox = boxLims(:, 2) - boxLims(:, 1); 
   for k = 1:numel(schedule.control)
      grd{k} = (dBox / objScaling) .* grd{k}; 
   end
end
    

