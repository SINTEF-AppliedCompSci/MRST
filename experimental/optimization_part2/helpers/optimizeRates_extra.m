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
   
   % Convergence details:
   opt.lineSearchMaxIt = 10;
   opt.gradTol         = 1e-4;
   opt.objChangeTol    = 1e-4;  
   
   % Penalization type:
   opt.penalize_type = 'leakage'; % 'leakage','leakage_at_infinity','pressure'
   
   % for default obj fun (leak_penalizer):
   opt.leak_penalty = 10; % Only used if 'obj_fun' not provided 
   
   % for pressure_penalizer:
   opt.pressure_penalty = [];
   opt.pressure_limit   = [];
   
   % for using leak_penalizer_at_infinity:
   %opt.penalize_leakage_at_infinity = false; % if true, uses alternative obj_fun to the default one
   opt.surface_pressure = [];
   opt.rho_water = [];
   
   % for accounting for well cost:
   opt.account_well_cost = false;
   opt.well_initial_cost = [];      % USD
   opt.well_operation_cost = [];    % USD/tonne
   opt.co2_tax_credit = [];         % USD/tonne
   opt.nonlinear_well_cost = true;    % if true, opt.alpha must not be empty
   opt.alpha = [];
   
   opt.last_control_is_migration = false; % if true, constrain last control
                                          % to zero rate
   opt.obj_fun = @(dummy) 0;
   opt = merge_options(opt, varargin{:});
   
   if opt.nonlinear_well_cost, assert(~isempty(opt.alpha)), end

   % default objective function
   opt.obj_fun = @(wellSols, states, schedule, varargin) ...
                  leak_penalizer(model, wellSols, states, schedule, opt.leak_penalty, ...
                                 varargin{:});
   
   % additional part of objective function (if opt.pressure_penalty
   % provided, obj fun will also penalize pressure close to pressure limit)
   if strcmpi(opt.penalize_type,'pressure')
       obj_funA = opt.obj_fun; clear opt.obj_fun; opt.obj_fun = @(dummy) 0;
       
       obj_funB = @(states, varargin) ...
          pressure_penalizer(model, states, schedule, opt.pressure_penalty, ...
                                    opt.pressure_limit, varargin{:});
       
       opt.obj_fun = @(wellSols, states, schedule, varargin) ...
           cell_subtract(obj_funA(wellSols, states, schedule, varargin{:}), ...
                         obj_funB(states, varargin{:}));                    
   end
   
   ta = trapAnalysis(model.G,false); % potentially used in evaluate_obj for spill-point dynamics
   
   % another possible objection function to use:
   % penalizes leakage which is determined at infinity, by computing the
   % volume of co2 accumulated (i.e., remaining) in formation at time
   % infinity using spill-point dynamics.
   % @@ requires testing...
   if strcmpi(opt.penalize_type,'leakage_at_infinity')
       assert(~isempty(opt.surface_pressure))
       assert(~isempty(opt.rho_water))
       opt.obj_fun = @(wellSols, states, schedule, varargin) ...
                  leak_penalizer_at_infinity(model, wellSols, states, ...
                  schedule, opt.leak_penalty, opt.surface_pressure, ...
                  opt.rho_water, ta, varargin{:});
   end
   
   % another possible term to add to objective function:
   % a term to account for well cost/operation
   if opt.account_well_cost
       assert(~isempty(opt.well_initial_cost));
       assert(~isempty(opt.well_operation_cost));
       assert(~isempty(opt.co2_tax_credit));
       assert(~isempty(opt.alpha));
       
       obj_funC = opt.obj_fun; clear opt.obj_fun; opt.obj_fun = @(dummy) 0;
       
       obj_funD = @(wellSols, schedule, varargin) account_well_cost(model, wellSols, schedule, ...
           opt.well_initial_cost, opt.well_operation_cost, opt.co2_tax_credit, opt.alpha, opt.nonlinear_well_cost, ...
           varargin{:});
       
%        opt.obj_fun = @(wellSols, states, schedule, varargin) ...
%            cell_subtract( cell_scalar_multiply( obj_funC(wellSols, states, schedule, varargin{:}), opt.co2_tax_credit*1e9 ), ...
%                          obj_funD(wellSols, schedule, varargin{:}) );
                     
       opt.obj_fun = @(wellSols, states, schedule, varargin) ...
           report_savings(  obj_funC(wellSols, states, schedule, varargin{:}), ...
                            opt.co2_tax_credit*1e9, ...
                            obj_funD(wellSols, schedule, varargin{:}), ...
                            varargin{:} );
   end
   
   
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
      %%
    figure(19),clf
    %plotCellData(model.G,init.states{end}.s(:,2)),colorbar;%NBhmn
    [h, h_max] = computePlumeHeight(model.G, init.states{end}, model.fluid.res_water,model.fluid.res_gas );
    %{
    G=model.G;
    G=mrstGridWithFullMappings(G);
    pcells=diff(G.n
    hn=accumarray(nocells,G.nodes.cells);
    nn=accumarray(nocells,G.nodes.cells
    %}
    plotCellData(model.G,h),colorbar;%NBhmn
    %%
      if strcmpi(opt.penalize_type,'pressure')
          init.obj_val_steps_A = cell2mat( obj_funA(init.wellSols, init.states, init.schedule) );
          init.obj_val_steps_B = cell2mat( obj_funB(init.states) );
          init.obj_val_steps = init.obj_val_steps_A - init.obj_val_steps_B;

          if ~ishandle(50)
            figure(50)
          else
            % Avoid stealing focus if figure already exists
            set(0, 'CurrentFigure', 50); clf(50)
          end
          hold on;
          plot(1:numel(init.obj_val_steps), init.obj_val_steps_A, 'x')
          plot(1:numel(init.obj_val_steps), -init.obj_val_steps_B, 'o')
          plot(1:numel(init.obj_val_steps), init.obj_val_steps, '+')
          legend('objA','-objB','objA - objB', 'Location','SouthWest')
          drawnow
          
      else
          init.obj_val_steps = cell2mat(opt.obj_fun(init.wellSols, ...
                                                init.states, ...
                                                init.schedule));
      end
      init.obj_val_total = sum(init.obj_val_steps);

      
      % Use value of objective function before optimization as scaling.
      % Otherwise, the scale of obj_val_steps_B can be much greater than _A
      opt.obj_scaling = abs(init.obj_val_total);
      if strcmpi(opt.penalize_type,'pressure')
        opt.obj_scaling = abs( sum(init.obj_val_steps_A) );
        
      %elseif opt.account_well_cost % @@ special scaling required ?
      %  opt.obj_scaling = abs( sum( cell2mat( obj_funC(init.wellSols, init.states, init.schedule) ) ) );
      
      end
   end
   
%     init.obj_val_steps_A_full = cell2mat(leak_penalizer_Rerun(model, init.wellSols, ...
%         init.states, schedule, opt.leak_penalty));
   
   
   %% Return if optimization should be skipped
   if opt.dryrun
      optim = init;
      history = [];
      return;
   end
   
   
   %% Define limits, scaling and objective function
   if strcmpi(schedule.control(1).W(1).type, schedule.control(2).W(1).type)
       % injection well-type is the same as migration well-type
       % thus same scalings apply to both control periods
       scaling.boxLims = [  min_wvals(:),     max_wvals(:)  ];
   else
       % injection well-type is different from migration well-type
       % thus use 2 separate scalings for control periods
       scaling.boxLims = [  min_wvals(:)    , max_wvals(:)    ; ...
                            min_mig_rates(:), min_mig_rates(:) + sqrt(eps)];
   end
   scaling.obj = opt.obj_scaling;
   
   obj_evaluator = @(u) evaluate_objective(u, opt.obj_fun, model, initState, ...
                                           schedule, scaling, opt.leak_penalty, ta); 

   %% Define constraints
   
   linEqS = [];
   if opt.last_control_is_migration
      % Constrain rates of last step to zero
      linEq = struct('A', eye(num_wells), 'b', min_mig_rates(:));
      
      if strcmpi(schedule.control(1).W(1).type, schedule.control(2).W(1).type)
        linEqS = setupConstraints(linEq, schedule, scaling);
        % keep only the relations pertaining to the last control step
        last_step_ix = num_wells * (numel(schedule.control) - 1) + 1;
        linEqS.A = linEqS.A(last_step_ix:end, :);
        linEqS.b = linEqS.b(last_step_ix:end);
      else
        % relations pertaining to last control step
        s2.control = schedule.control(2);
        scaling2.boxLims = scaling.boxLims(end/2+1:end,:);
        linEqS = setupConstraints(linEq, s2, scaling2);
        linEqS.A = [zeros(size(linEqS.A)), linEqS.A];
      end
   end
   
   
   %% Call optimization routine
   if strcmpi(schedule.control(1).W(1).type, schedule.control(2).W(1).type)
        u = schedule2control(schedule, scaling);
   else
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
   end
   
   [~, u_opt, history] = unitBoxBFGS(u, obj_evaluator, 'linEq', linEqS, ...
       'lineSearchMaxIt', opt.lineSearchMaxIt, ...
       'gradTol'        , opt.gradTol, ...
       'objChangeTol'   , opt.objChangeTol,...
       'useBFGS',true);
   history.linEqS=linEqS;
   %% Preparing solution structures
   if strcmpi(schedule.control(1).W(1).type, schedule.control(2).W(1).type)
       optim.schedule = control2schedule(u_opt, schedule, scaling);
   else
       % copy schedule of control.BC and step
       optim.schedule = schedule;
       % get optimized schedule (of Wells) from optimized control vector
       s1 = control2schedule(u_opt(1:end/2), s1, scaling1);
       s2 = control2schedule(u_opt(end/2+1:end), s2, scaling2);
       optim.schedule.control(1) = s1.control;
       optim.schedule.control(2) = s2.control;
   end
   
   [optim.wellSols, optim.states] = simulateScheduleAD(initState, ...
                                                       model, ...
                                                       optim.schedule);
   
   optim.obj_val_steps = opt.obj_fun(optim.wellSols, optim.states, optim.schedule);
   optim.obj_val_total = sum(cell2mat(optim.obj_val_steps));
end

% ----------------------------------------------------------------------------

function c = cell_subtract(a, b)
    % we assume that a and b should be single-indexed
    c = a;
    for i = 1:numel(c)
        c{i} = c{i} - b{i};
    end
end
 
function c = cell_scalar_multiply(a, scalar)
    c = a;
    for i = 1:numel(c)
        c{i} = c{i}*scalar; 
    end
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
      sF = state.s(:,1);   % for dissolution
      if isfield(model.fluid,'rsSat') % dissolution is true
        sGmax = state.sGmax; % for dissolution
        rs = state.rs;       % for dissolution
      end
      if opt.ComputePartials
         if isfield(model.fluid,'rsSat') % dissolution is true
            [p, sG, sGmax, rs, qWs, qGs, pBHP] = initVariablesADI(p, sG, sGmax, rs, qWs, qGs, pBHP); 
         else
            %[p, sG, pBHP, qWs, qGs] = initVariablesADI(p, sG, pBHP, qWs, qGs);%#ok
            [p, sG, qWs, qGs, pBHP] = initVariablesADI(p, sG, qWs, qGs, pBHP);%#ok
         end
      end
      dt = dts(step);
      injInx = (vertcat(sol.sign) > 0);
      obj{step} = dt * spones(ones(1, nW)) * ((1-penalty) * injInx .* qGs);

      krull = krull +  dt * spones(ones(1, nW)) * ( injInx .* qGs);
      
      if (tSteps(step) == num_timesteps)
         bG = model.fluid.bG(p);
         bW = model.fluid.bW(p);
         if ~isfield(model.rock,'ntg')
             model.rock.ntg = ones(model.G.cells.num,1); % in case ntg doesn't exist
         end 
         pvol = model.G.cells.volumes .* model.G.cells.H .* model.rock.poro .* model.rock.ntg;
         vol = ones(1, model.G.cells.num) * (pvol .* model.fluid.pvMultR(p) .* bG .* sG);
         if isfield(model.fluid,'rsSat') % dissolution is true
             vol_diss = ones(1, model.G.cells.num) * (pvol .* model.fluid.pvMultR(p) .* rs .* bW .* sF);
             vol = vol + vol_diss;
         end
         obj{step} = obj{step} + penalty * vol;
         
         if ~opt.ComputePartials
            fprintf('Total injected: %f (m3)\n', double(krull));
            fprintf('Total leaked: %f (m3)\n', double(krull - vol));
            fprintf('Penalty: %f \n', penalty)
            fprintf('Score: %f (m3)\n\n', double(krull) - penalty * double(krull-vol));
         end
      end
      obj{step} = obj{step} * model.fluid.rhoGS/1e12; % converted to Gt
   end
end

% ----------------------------------------------------------------------------

function obj = leak_penalizer_at_infinity(model, wellSols, states, schedule, ...
    penalty, surf_press, rho_water, ta, varargin)
% new objective function to penalize leakage 'felt' at infinity, using
% spill-point dynamics.

% pressure is assumed to be hydrostatic at time infinity.

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
      p_future = compute_hydrostatic_pressure(model.G, rho_water, surf_press);
      sG = state.s(:,2);
      if opt.ComputePartials
         %[p, sG, pBHP, qWs, qGs] = initVariablesADI(p, sG, pBHP, qWs, qGs);%#ok
         [p, sG, sGmax, qWs, qGs, pBHP] = initVariablesADI(p, sG, sGmax, qWs, qGs, pBHP);%#ok
         %state0 = states{tSteps(step-1)}
         %sGmax0 = state0.sGmax;
         %sGmax = max(sG,sGmax0);
      end
      dt = dts(step);
      injInx = (vertcat(sol.sign) > 0);
      % calculate the "(1 - C) M^inj" part of obj fun J
      obj{step} = dt * spones(ones(1, nW)) * ((1-penalty) * injInx .* qGs);

      krull = krull +  dt * spones(ones(1, nW)) * ( injInx .* qGs);
      
      if (tSteps(step) == num_timesteps)
         bG = model.fluid.bG(p);
         if ~isfield(model.rock,'ntg')
             model.rock.ntg = ones(model.G.cells.num,1); % in case ntg doesn't exist
         end 
         pvol = model.G.cells.volumes .* model.G.cells.H .* model.rock.poro .* model.rock.ntg;
         vol = ones(1, model.G.cells.num) * (pvol .* model.fluid.pvMultR(p) .* bG .* sG);
         
         % calculate the "C M^a" part of obj fun J, where M^a is the mass
         % accumulated (or remaining) at t = infinity. To determine this
         % mass, we calculate the vol expected to be retained using
         % spill-pt dynamics.
%          vol_inf = vol_at_infinity( model.G, ...
%              model.rock.poro .* model.fluid.pvMultR(state.pressure), ...
%              state.s(:,2) .* model.fluid.bG(state.pressure), ...
%              model.fluid.res_water );
%          rock2D.ntg  = ones(model.G.cells.num,1);
%          if isfield(model.rock,'ntg')
%              rock2D.ntg = model.rock.ntg; % in case ntg doesn't exist
%          end 
%          rock2D.poro = model.rock.poro .* model.fluid.pvMultR(p);
%          vol_inf = vol_at_infinity( model.G, ...
%              rock2D, ... % or rock2D.poro .* model.fluid.pvMultR(p) ? @@
%              sG .* model.fluid.bG(p), ...
%              model.fluid.res_water, model.fluid.res_gas); % using possible ADI variables
         vol_inf = vol_at_infinity( model.G, ...
             model.rock, sG, model.fluid, p, p_future, 'ta', ta );
         % NB: vol_inf is at ref. depth
         % NB: sG, p may be ADI variables, thus vol_inf may also be ADI
         
         % J = (1-C) M^i + C M^a
         obj{step} = obj{step} + penalty * vol_inf;
         if ~opt.ComputePartials
            fprintf('Total injected: %f (m3)\n', double(krull));
            fprintf('Total leaked (by infinity): %f (m3)\n', double(krull - vol_inf));
            fprintf('Penalty: %f \n', penalty)
            fprintf('Score: %f (m3)\n\n', double(krull) - penalty * double(krull-vol_inf));
         end
      end
      obj{step} = obj{step} * model.fluid.rhoGS/1e12; % vol * rho, in Gt.
   end


end

% ----------------------------------------------------------------------------

function obj = account_well_cost(model, wellSols, schedule, initial_cost, operation_rate, tax_credit_rate, alpha, nonlinear_well_cost, varargin )
% operation_rate is in units of USD per tonne CO2
% tax_credit_rate is in units of USD per tonne CO2
% initial_cost is in units of USD
% obj{steps} returned is in units of USD

   opt.ComputePartials = false;
   opt.tStep = [];
   opt = merge_options(opt, varargin{:});
   
   num_timesteps = numel(schedule.step.val);
   num_inj_timesteps = numel(schedule.step.val(schedule.step.control==1));
   tSteps = opt.tStep;
   if isempty(tSteps)
      numSteps = numel(wellSols); % @@ use wellSols or schedule?
      tSteps = (1:numSteps)';
      dts = schedule.step.val;
   else
      assert(numel(tSteps) == 1);
      numSteps = 1;
      dts = schedule.step.val(opt.tStep);
   end
   
   obj = repmat({[]}, numSteps, 1);
   total_inj = zeros(numSteps, numel(wellSols{1,1}));
   obj_per_well = zeros(numSteps, numel(wellSols{1,1}));

   for step = 1:numSteps
      sol = wellSols{tSteps(step)};
      nW = numel(sol);
      pBHP = zeros(nW, 1); % place holder
      qGs = vertcat(sol.qGs);
      qWs = vertcat(sol.qWs);
      p = zeros(model.G.cells.num,1); % place holder
      sG = zeros(model.G.cells.num,1); % place holder
      if opt.ComputePartials
         %[p, sG, pBHP, qWs, qGs] = initVariablesADI(p, sG, pBHP, qWs, qGs);%#ok
         [p, sG, qWs, qGs, pBHP] = initVariablesADI(p, sG, qWs, qGs, pBHP);%#ok
      end
      dt = dts(step);
      injInx = (vertcat(sol.sign) > 0);
      
      % For each time step, compute the well cost
      
      %cost_for_having = ones(nW, 1) * (initial_cost/num_inj_timesteps); % USD
      cost_for_having = ones(nW, 1) * (initial_cost/num_inj_timesteps) .* max(0, sign(qGs - sqrt(eps))); % @@ assuming sqrt(eps) is min val
                    
      % NB: operation_cost is [USD/tonne] (1 tonne = 1e3 kg)
      % NB: dt*qGs is [m3] at ref depth
      % NB: dt*qGs*rhoGS is [kg]
      cost_for_operating = operation_rate * dt * injInx .* qGs .* model.fluid.rhoGS/1e3 .* max(0, sign(qGs - sqrt(eps)));     % USD
                       
      
      if ~nonlinear_well_cost
          
        % Linear total well cost:
        obj{step} = sum( (cost_for_having + cost_for_operating) .* max(0, sign(qGs - 2*sqrt(eps))) ); % USD
        % for plotting:
        obj_per_well(step,1:nW) = cost_for_having + cost_for_operating; % USD
      
      else
          
        % Non-linear total well cost:
        % Critical injection mass rate (for this time step) to make well's investment cost worthwhile:
        %q_crit = (initial_cost/num_inj_timesteps)/(operation_rate) * 1e3 / model.fluid.rhoGS / dt; % m3/s
        q_crit = (initial_cost/num_inj_timesteps)/(tax_credit_rate - operation_rate) * 1e3 / model.fluid.rhoGS / dt; % m3/s
        obj{step} = sum( cost_for_having .* tanh(qGs./(alpha*q_crit)) + cost_for_operating ); % USD
        % for plotting:
        obj_per_well(step,1:nW) = cost_for_having .* tanh(qGs./(alpha*q_crit)) + cost_for_operating; % USD
      
      end
      
      % for plotting purposes
      total_inj(step,1:nW) = qGs .* dt .* model.fluid.rhoGS .* max(0, sign(qGs - sqrt(eps))); % kg

      
      if (tSteps(step) == num_timesteps)
         if ~opt.ComputePartials
            fprintf('Well cost is %f USD. \n\n', sum(cell2mat(obj)) );
            %
            % Total well cost
            if ~ishandle(12)
                figure(12)
            else
                % Avoid stealing focus if figure already exists
                set(0, 'CurrentFigure', 12); clf(12);
            end
            plot_total_well_cost_function( 12, initial_cost, operation_rate, tax_credit_rate, alpha ) % opens new figure
            hold on; plot(sum(total_inj)/1e9, sum(obj_per_well)/1e6, 'xk')
            text(sum(total_inj)/1e9, sum(obj_per_well)/1e6, cellstr(num2str([1:nW]')))
            %
            % Max savings possible per well (assuming no leakage)
            if ~ishandle(13)
                figure(13)
            else
                % Avoid stealing focus if figure already exists
                set(0, 'CurrentFigure', 13); clf(13);
            end
            plot_total_savings_function( 13, initial_cost, operation_rate, tax_credit_rate, alpha )
            hold on; plot(sum(total_inj)/1e9, tax_credit_rate .* sum(total_inj)/1e9 - sum(obj_per_well)/1e6, 'xk')
            text(sum(total_inj)/1e9, tax_credit_rate .* sum(total_inj)/1e9 - sum(obj_per_well)/1e6, cellstr(num2str([1:nW]')))
         end
      end
   end

end


% ----------------------------------------------------------------------------

function obj = pressure_penalizer(model, states, schedule, penalty, plim, varargin)
% states.pressure is a cell array.
% schedule is only used for time steps.
% penalty is a scalar.
% plim is a cell array.


% format of objective function:
%   obj = max(0, sign(p - plim)) * penalty * (p - plim)^2
% obj is computed for each time step (for both inj and mig periods).


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
   k   = 2;
   max_amount_surp      = zeros(numSteps, 1);
   max_amount_surp_cinx = zeros(numSteps, 1);
   min_amount_under     = zeros(numSteps, 1);
   min_amount_under_cinx= zeros(numSteps, 1);
   for step = 1:numSteps     
      state = states{tSteps(step)}; %@@ +1?      
      p = state.pressure;
      % keep track of amount over or amount under plim at each time step
      max_amount_surp(step) = max(0,max(p-plim));
      if max_amount_surp(step) > 0
          [~,cinx] = max(p-plim);
          max_amount_surp_cinx(step) = cinx;
      end
      min_amount_under(step) = max(0,min(plim-p));
      if min_amount_under(step) > 0
         [~,cinx] = min(max(0,(plim-p)));
         min_amount_under_cinx(step) = cinx;
      end
      if opt.ComputePartials
        sG = state.s(:,2);   % place holders
        nW = numel(schedule.control(1).W);
        pBHP = zeros(nW, 1); % place holders
        qGs = pBHP;          % place holders
        qWs = pBHP;          % place holders
        [p, ~, ~, ~, ~] = initVariablesADI(p, sG, qWs, qGs, pBHP); 
      end
      dt = dts(step);
      tmp = max(0, sign(p - plim)) .* penalty .* (p - plim).^k/1e12; % @@ scaling to MPa
      tmp = tmp .* model.G.cells.volumes;
      obj{step} = sum( tmp )./sum(model.G.cells.volumes);
      obj{step} = obj{step} * dt; % @@ if uncommented, cp values should be very small
      if (tSteps(step) == num_timesteps)
      % no need to compute another portion of the obj fun here
         if ~opt.ComputePartials
             msg1 = 0;
             msg2 = 0;
             if any(max_amount_surp > 0)
                 [val,tinx] = max(max_amount_surp);
                 cinx = max_amount_surp_cinx(tinx);
                 msg1 = 100*val/plim(cinx);
             else
                 [val,tinx] = min(min_amount_under);
                 cinx = min_amount_under_cinx(tinx);
                 msg2 = 100*val/plim(cinx);
             end
             fprintf('Surpassed Plimit by %f (percent) of Plimit.\n', msg1)
             fprintf('Approached Plimit by %f (percent) of Plimit.\n', msg2)
         end
      end

   end
end

% ----------------------------------------------------------------------------

function [val, der, wellSols, states] = ...
       evaluate_objective(u, obj_fun, model, initState, schedule, scaling, cp, ta) 
% cp and ta passed in for plotting purposes
   
   if ~strcmpi(schedule.control(1).W(1).type, schedule.control(2).W(1).type)
       % assume:
       % first half of controls are for injection period (could be bhp or rate)
       % second half of controls are for migration period (should be rate)
       s1.control = schedule.control(1);
       s2.control = schedule.control(2);
       scaling1.boxLims = scaling.boxLims(1:end/2,:);
       scaling2.boxLims = scaling.boxLims(end/2+1:end,:);
       % where u = [u1; u2]
   end

   minu = min(u);
   maxu = max(u);
   if or(minu < -eps , maxu > 1+eps)
      warning('Controls are expected to lie in [0 1]\n')
   end
   
   if strcmpi(schedule.control(1).W(1).type, schedule.control(2).W(1).type)
       boxLims = scaling.boxLims;
   else
       boxLims1 = scaling1.boxLims;
       boxLims2 = scaling2.boxLims;
   end
   if isfield(scaling, 'obj')
      objScaling = abs(scaling.obj);
   else
      objScaling = 1;
   end
   
   % update schedule:
   if strcmpi(schedule.control(1).W(1).type, schedule.control(2).W(1).type)
       schedule = control2schedule(u, schedule, scaling);
   else
       s1 = control2schedule(u(1:end/2), s1, scaling1);
       s2 = control2schedule(u(end/2+1:end), s2, scaling2);
       schedule.control(1) = s1.control;
       schedule.control(2) = s2.control;
   end
   
   % run simulation:
   [wellSols, states] = simulateScheduleAD(initState, model, schedule);
   %figure(20),clf
   figure(),clf
   plotCellData(model.G,states{end}.s(:,2)),colorbar;%NBhmn
   % compute objective:
   vals = obj_fun(wellSols, states, schedule); % @@ ensure wellSols{}.sign is correct
   if ~ishandle(11)
     figure(11)
   else
     % Avoid stealing focus if figure already exists
     set(0, 'CurrentFigure', 11); clf(11);
   end
   h1 = plot(1:numel(vals), [vals{:}], '+');
   val = sum(cell2mat(vals))/abs(objScaling);
   %legend([h1],{['scaled obj val: ',num2str(val)]},'Location','SouthWest')
   title(['scaled obj val: ',num2str(val)])
   drawnow

   % visualize:
   % assuming surface pressure = 1 * atm
   %plotObjValues(model, wellSols, states, schedule, cp, 1 * atm, ta)
   
   % run adjoint:
   if nargout > 1
      objh = @(tstep)obj_fun(wellSols, states, schedule, 'ComputePartials', true, 'tStep', tstep);
      g    = computeGradientAdjointAD(initState, states, model, schedule, objh);
      % scale gradient:
      if strcmpi(schedule.control(1).W(1).type, schedule.control(2).W(1).type)
          der = scaleGradient(g, schedule, boxLims, objScaling);
          der = vertcat(der{:});
      else
          der1 = scaleGradient({g{1}}, s1, boxLims1, objScaling);
          der2 = scaleGradient({g{2}}, s2, boxLims2, objScaling);
          der = vertcat(der1{:}, der2{:});
      end
%       %% @@ 
%       % Compute numeric derivative, to verify gradient
%       vd = u*0;
%       du = 1e-7;
%       for i = 1:numel(u)
%          u_tmp = u;
%          u_tmp(i) = u_tmp(i) + du; % to compute partial derivative along i
%          tmp_schedule = control2schedule(u_tmp, schedule, scaling);
%          [ws, st] = simulateScheduleAD(initState, model, tmp_schedule);
%          tmp_val = obj_fun(ws, st, tmp_schedule);
%          tmp_val = sum(cell2mat(tmp_val))/abs(objScaling);
%          vd(i) = (tmp_val - val) / du;
%       end
%       vd;
   end
end

function obj = report_savings(  J1, co2_tax_credit, J2, varargin )
% J1: obj values [Gt] for each time step, computed by an obj function
% (either 'leak_penalizer', 'leak_penalizer_at_infinity', 'pressure_penalizer')
% co2 tax credit [USD/Gt]
% J2: obj values [USD] for each time step, computed by account_well_cost

   opt.ComputePartials = false;
   opt.tStep = [];
   opt = merge_options(opt, varargin{:});

    % simply compute amount USD saved (neg implies USD spent)
    % obj = MassStored * tax_credit - sum_i ( InitialCost_i + OperationCost_i ); for i=1:num_wells
    obj = cell_subtract( cell_scalar_multiply(J1,co2_tax_credit), J2 );
    
    if ~opt.ComputePartials
        % reporting:
        fprintf('\nMass stored (after penalization): %f Gt\n', sum(full(cell2mat(J1))) ) % @@ J1 may be sparse if cp=1
        fprintf('Value of this mass stored: %f USD\n', sum( full(cell2mat(cell_scalar_multiply(J1,co2_tax_credit))) ) )
        fprintf('Well costs: %f USD\n', sum(cell2mat(J2)) )  
        fprintf('Savings: %f USD\n\n', sum(cell2mat(obj)) )
    end
end


% ----------------------------------------------------------------------------

function grd = scaleGradient(grd, schedule, boxLims, objScaling)
   dBox = boxLims(:, 2) - boxLims(:, 1); 
   for k = 1:numel(schedule.control)
      grd{k} = (dBox / objScaling) .* grd{k}; 
   end
end

% ----------------------------------------------------------------------------

function p = compute_hydrostatic_pressure(Gt, rho_water, surface_pressure)

    p = rho_water * norm(gravity()) * Gt.cells.z + surface_pressure;
end

function plotObjValues(model, wellSols, states, schedule, cp, surf_press, ta)

    % Obtain the obj values: (Xyrs is the number of migration years)
    % Also, mass or volume inventories can be shown
    [obj_val_steps_Xyrs, ~, Mi_tot, Ma] = leak_penalizer_Rerun(model, wellSols, ...
                        states, schedule, cp); 
    [obj_val_steps_Xyrs_noPenalty] = leak_penalizer_Rerun(model, wellSols, ...
                        states, schedule, 1);

    [obj_val_steps_Xyrs_future, Mi_tot_inf, Ma_inf] = leak_penalizer_at_infinity_Rerun(model, ...
                        wellSols, states, schedule, cp, surf_press, model.fluid.rhoWS, ta);
    [obj_val_steps_Xyrs_future_noPenalty, ~, ~] = leak_penalizer_at_infinity_Rerun(model, ...
                        wellSols, states, schedule, 1, surf_press, model.fluid.rhoWS, ta);
    % NB: Mi_tot and Mi_tot_inf should be the same given the same injection schedule 
    assert(all(Mi_tot == Mi_tot_inf))
    
    % Compare obj values:
    if ~ishandle(24)
     figure(24)
    else
     % Avoid stealing focus if figure already exists
     set(0, 'CurrentFigure', 24); clf(24);
    end
    hold on; set(gcf,'Position',[2749 166 1120 411])
    plot(convertTo(cumsum(schedule.step.val), year), ...
        [obj_val_steps_Xyrs{:}], 'o') % Gt
    plot(convertTo(cumsum(schedule.step.val), year), ...
        [obj_val_steps_Xyrs_future{:}], 'o')
    hl = legend('leakage penalized','future leakage penalized');
    set(hl,'Location','best')
    ylabel('J(t), (Gt CO2)'); xlabel('time (years)')
    title(['C=',num2str(cp)])
    box; grid;
    set(gca,'FontSize',16)
    
    % Compare obj values with c=1 (which corresponds to mass inventory):
    if ~ishandle(25)
     figure(25)
    else
     % Avoid stealing focus if figure already exists
     set(0, 'CurrentFigure', 25); clf(25);
    end
    hold on; set(gcf,'Position',[3873 165 1119 412])
    plot(convertTo(cumsum(schedule.step.val), year), ...
        [obj_val_steps_Xyrs_noPenalty{:}], 'o') % Gt
    plot(convertTo(cumsum(schedule.step.val), year), ...
        [obj_val_steps_Xyrs_future_noPenalty{:}], 'o')
    plot(convertTo(cumsum(schedule.step.val), year), Mi_tot , 'x')
    plot(convertTo(cumsum(schedule.step.val), year), Ma , '+')
    plot(convertTo(cumsum(schedule.step.val), year), Ma_inf , '+')
    hl = legend('leakage penalized','future leakage penalized', ...
        'tot Mi(t) (mass injected)', 'Ma(t) (mass remaining)', 'future Ma(t) (future mass remaining)');
    set(hl,'Location','best')
    ylabel('J(t), (Gt CO2)'); xlabel('time (years)')
    title('C=1')
    box; grid;
    set(gca,'FontSize',16)

end
    

