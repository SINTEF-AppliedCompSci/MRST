function [optim, init, history] = optimizeRates(initState, model, schedule, ...
                                                min_rates, max_rates, varargin)
   
% optim.obj_val_steps
% optim.obj_val_total
% optim.schedule
% optim.wellSols
% optim.states

% opt.obj_fun should have the following interface
% function obj = obj_fun(wellSols, states, schedule, varargin);
% Where 'varargin' should support 'ComputePartials' and 'tStep'
% 'obj' should here be a cell array with objective values per timestep.  If
% 'ComputePartials' is given, then the values should be ADI.
% @@ Do ADI variables need to be exactly as in the equation??


   moduleCheck('ad-core', 'ad-props', 'optimization');

   opt.obj_scaling = [];  % Compute explicitly if not provided from outside
   opt.leak_penalty = 10; % Only used if 'obj_fun' not provided 
   opt.last_control_is_migration = false; % if true, constrain last control to zero rate

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
      opt.obj_scaling = abs(init.obj_val_total);
   end
   
   %% Define limits, scaling and objective function
   
   scaling.boxLims = [min_rates(:), max_rates(:)];
   scaling.obj     = opt.obj_scaling;
   
   obj_evaluator = @(u) evaluate_objective(u, opt.obj_fun, model, initState, ...
                                           schedule, scaling); 

   %% Define constraints
   
   linEqS = [];
   if opt.last_control_is_migration
      % Constrain rates of last step to zero
      linEq = struct('A', eye(num_wells), 'b', min_rates(:));
      linEqS = setupConstraints(linEq, schedule, scaling);
      
      % keep only the relations pertaining to the last control step
      last_step_ix = num_wells * (numel(schedule.control) - 1) + 1;
      linEqS.A = linEqS.A(last_step_ix:end, :);
      linEqS.b = linEqS.b(last_step_ix:end);
   end
   
   
   %% Call optimization routine
   
   u = schedule2control(schedule, scaling);
   [~, u_opt, history] = unitBoxBFGS(u, obj_evaluator, 'linEq', linEqS, 'lineSearchMaxIt', 20);
   
   %% Preparing solution structures
   
   optim.schedule = control2schedule(u_opt, schedule, scaling);

   [optim.wellSols, optim.states] = simulateScheduleAD(initState, ...
                                                       model, ...
                                                       optim.schedule);
   
   optim.obj_val_steps = opt.obj_fun(optim.wellSols, optim.states, optim.schedule);
   optim.obj_val_total = sum(cell2mat(optim.obj_val_steps));
end

% ----------------------------------------------------------------------------

function obj = leak_penalizer(model, wellSols, states, schedule, penalty, varargin)
   
   opt.ComputePartials = false;
   opt.tStep = [];
   opt = merge_options(opt, varargin{:});
   
   num_timesteps = numel(schedule.step.val);
   tSteps = opt.tStep;
   if isempty(tSteps)
      numSteps = numel(states);
      tSteps = (1:numSteps)';
      dts = schedule.step.val;
   else
      assert(numel(tSteps) == 1);
      numSteps = 1;
      dts = schedule.step.val(opt.tStep);
   end
   
   obj = repmat({[]}, numSteps, 1);
   krull = 0; % @@
   for step = 1:numSteps
      sol = wellSols{tSteps(step)};
      state = states{tSteps(step)}; %@@ +1?
      nW = numel(sol);
      pBHP = zeros(nW, 1); % place holder
      qGs = vertcat(sol.qGs);
      qWs = vertcat(sol.qWs);
      p = state.pressure;
      sG = state.s(:,2);
      if opt.ComputePartials
         [p, sG, pBHP, qWs, qGs] = initVariablesADI(p, sG, pBHP, qWs, qGs);%#ok
      end
      dt = dts(step);
      injInx = (vertcat(sol.sign) > 0);
      obj{step} = dt * spones(ones(1, nW)) * ((1-penalty) * injInx .* qGs);

      krull = krull +  dt * spones(ones(1, nW)) * ( injInx .* qGs);
      
      if (tSteps(step) == num_timesteps)
         bG = model.fluid.bG(p);
         pvol = model.G.cells.volumes .* model.G.cells.H .* model.rock.poro;      
         vol = ones(1, model.G.cells.num) * (pvol .* model.fluid.pvMultR(p) .* bG .* sG);
         obj{step} = obj{step} + penalty * vol;
         if ~opt.ComputePartials
            fprintf('Total injected: %f\n', double(krull));
            fprintf('Total leaked: %f\n', double(krull - vol));
            fprintf('Score: %f\n\n', double(krull) - penalty * double(krull-vol));
         end
      end
      obj{step} = obj{step} * model.fluid.rhoGS/1e12; %@@
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
      objScaling = abs(scaling.obj);
   else
      objScaling = 1;
   end
   
   % update schedule:
   schedule = control2schedule(u, schedule, scaling);
   
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
      der = vertcat(der{:});
      
      % %% @@ 
      % % Compute numeric derivative, to verify gradient
      % vd = u*0;
      % du = 1e-7;
      % for i = 1:numel(u)
      %    u_tmp = u;
      %    u_tmp(i) = u_tmp(i) + du; % to compute partial derivative along i
      %    tmp_schedule = control2schedule(u_tmp, schedule, scaling);
      %    [ws, st] = simulateScheduleAD(initState, model, tmp_schedule);
      %    tmp_val = obj_fun(ws, st, tmp_schedule);
      %    tmp_val = sum(cell2mat(tmp_val))/abs(objScaling);
      %    vd(i) = (tmp_val - val) / du;
      % end
      % vd;
   end
end

% ----------------------------------------------------------------------------

function grd = scaleGradient(grd, schedule, boxLims, objScaling)
   dBox = boxLims(:, 2) - boxLims(:, 1); 
   for k = 1:numel(schedule.control)
      grd{k} = (dBox / objScaling) .* grd{k}; 
   end
end
    

