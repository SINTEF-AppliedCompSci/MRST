function [optim, init, history] = optimizeRates_extra(initState, model, schedule, ...
                                                min_wvals, max_wvals, min_mig_rates, ...
                                                varargin)
   
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

   opt.dryrun = false;    % if 'true', no optimization will take place
   opt.obj_scaling = [];  % Compute explicitly if not provided from outside
   opt.leak_penalty = 10; % Only used if 'obj_fun' not provided 
   opt.last_control_is_migration = false; % if true, constrain last control
                                          % to zero rate
   opt.obj_fun = @(dummy) 0;
   opt = merge_options(opt, varargin{:});

   opt.obj_fun = @(wellSols, states, schedule, varargin) ...
                  leak_penalizer(model, wellSols, states, schedule, opt.leak_penalty, ...
                                 varargin{:});
   
   opt = merge_options(opt, varargin{:});
   
   num_wells = numel(schedule.control(1).W);
   assert(numel(max_wvals) == num_wells);
   assert(numel(min_wvals) == num_wells);
   assert(numel(min_mig_rates) == num_wells);
   
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
   
   %% Return if optimization should be skipped
   if opt.dryrun
      optim = init;
      history = [];
      return;
   end
   
   
   %% Define limits, scaling and objective function
   
   scaling.boxLims = [  min_wvals(:)    , max_wvals(:)    ; ...
                        min_mig_rates(:), min_mig_rates(:) + sqrt(eps)];
%    scaling1.boxLims = [min_wvals(:), max_wvals(:)];
%    scaling2.boxLims = [min_mig_rates(:), min_mig_rates(:) + sqrt(eps)];
%    scaling1.obj = opt.obj_scaling;
%    scaling2.obj = opt.obj_scaling;
   scaling.obj     = opt.obj_scaling;
   
   obj_evaluator = @(u) evaluate_objective(u, opt.obj_fun, model, initState, ...
                                           schedule, scaling); 

   %% Define constraints
   
   linEqS = [];
   if opt.last_control_is_migration
      % Constrain rates of last step to zero
      linEq = struct('A', eye(num_wells), 'b', min_mig_rates(:));
      %linEqS = setupConstraints(linEq, schedule, scaling);
      
      % keep only the relations pertaining to the last control step
%       last_step_ix = num_wells * (numel(schedule.control) - 1) + 1;
%       linEqS.A = linEqS.A(last_step_ix:end, :);
%       linEqS.b = linEqS.b(last_step_ix:end);
      
      % relations pertaining to last control step
      s2.control = schedule.control(2);
      scaling2.boxLims = scaling.boxLims(end/2+1:end,:);
      linEqS = setupConstraints(linEq, s2, scaling2);
      linEqS.A = [zeros(size(linEqS.A)), linEqS.A];
   end
   
   
   %% Call optimization routine
   
   %u = schedule2control(schedule, scaling);
   
   % Convert schedule-params (i.e., inj 'bhp', mig 'rate') to control
   % vector: 
   % - control(1) is bounded by bhp, control(2) is bounded by rates
   s1.control = schedule.control(1);
   s2.control = schedule.control(2);
   scaling1.boxLims = scaling.boxLims(1:end/2,:);
   scaling2.boxLims = scaling.boxLims(end/2+1:end,:);
   u1 = schedule2control(s1, scaling1);
   u2 = schedule2control(s2, scaling2);
   u = [u1; u2];
   
   % Convert schedule-params (i.e., inj 'bhp', mig 'rate') to control vector % @@
%    ui = [schedule.control(1).W.val]';
%    um = [schedule.control(2).W.val]';
%    u = [ui; um];
%    [umin, umax] = deal(scaling.boxLims(:,1), scaling.boxLims(:,2));
%    u  = (u-umin)./(umax-umin); 
   
   [~, u_opt, history] = unitBoxBFGS(u, obj_evaluator, 'linEq', linEqS, ...
       'lineSearchMaxIt', 20, 'gradTol',5e-3);
   
   %% Preparing solution structures
   
   %optim.schedule = control2schedule(u_opt, schedule, scaling);
   % copy schedule of control.BC and step
   optim.schedule = schedule;
   % get optimized schedule (of Wells) from optimized control vector
   s1 = control2schedule(u_opt(1:end/2), s1, scaling1);
   s2 = control2schedule(u_opt(end/2+1:end), s2, scaling2);
   optim.schedule.control(1) = s1.control;
   optim.schedule.control(2) = s2.control;
   
   

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
         %[p, sG, pBHP, qWs, qGs] = initVariablesADI(p, sG, pBHP, qWs, qGs);%#ok
         [p, sG, qWs, qGs, pBHP] = initVariablesADI(p, sG, qWs, qGs, pBHP);%#ok
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
   
   % assume:
   % first half of controls are for injection period (could be bhp or rate)
   % second half of controls are for migration period (should be rate)
    s1.control = schedule.control(1);
    s2.control = schedule.control(2);
    scaling1.boxLims = scaling.boxLims(1:end/2,:);
    scaling2.boxLims = scaling.boxLims(end/2+1:end,:);
    
    % where u = [u1; u2]
    
   minu = min(u);
   maxu = max(u);
   if or(minu < -eps , maxu > 1+eps)
      warning('Controls are expected to lie in [0 1]\n')
   end

%    boxLims = scaling.boxLims;
   boxLims1 = scaling1.boxLims;
   boxLims2 = scaling2.boxLims;
   if isfield(scaling, 'obj')
      objScaling = abs(scaling.obj);
   else
      objScaling = 1;
   end
   
   % update schedule:
   %schedule = control2schedule(u, schedule, scaling); %@@ use different scaling for bhp and rate controls
   
   s1 = control2schedule(u(1:end/2), s1, scaling1);
   s2 = control2schedule(u(end/2+1:end), s2, scaling2);
   schedule.control(1) = s1.control;
   schedule.control(2) = s2.control;
   
   % run simulation:
   [wellSols, states] = simulateScheduleAD(initState, model, schedule);
   
   % compute objective:
   vals = obj_fun(wellSols, states, schedule); % @@ ensure wellSols{}.sign is correct
   val  = sum(cell2mat(vals))/abs(objScaling);

   % run adjoint:
   if nargout > 1
      objh = @(tstep)obj_fun(wellSols, states, schedule, 'ComputePartials', true, 'tStep', tstep);
      g    = computeGradientAdjointAD(initState, states, model, schedule, objh);
      % scale gradient:
      %der = scaleGradient(g, schedule, boxLims, objScaling);
      der1 = scaleGradient({g{1}}, s1, boxLims1, objScaling);
      der2 = scaleGradient({g{2}}, s2, boxLims2, objScaling);
      %der = vertcat(der{:});
      der = vertcat(der1{:}, der2{:});
      
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
    

