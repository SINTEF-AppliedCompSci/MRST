function [optim, init, history] = optimizeControls(initState, model, schedule, ...
                                                min_wvals, max_wvals, varargin)
%
% Compute an optimal set of well controls ('rate' or 'bhp') for a proposed
% injection/migration scenario.
% 
% SYNOPSIS:
%   function [optim, init, history] = ...
%      optimizeControls(initState, model, schedule, min_wvals, max_wvals, varargin)
%
% DESCRIPTION:
%   The injection phase can consist of multiple control periods. If a
%   migration phase is included, it must be only one control period, and it
%   must be the last control period. Individual wells must be of the same
%   well type (either 'rate' or 'bhp') during the entire injection phase.
%   If there's a migration phase, all wells must be of well type 'rate'.
%   Individual well signs must remain fixed over all periods (+1 for
%   injector, -1 for producer).
%
%
% PARAMETERS:
%   initState - initial state
%   model     - simulation model (instance of CO2VEBlackOilTypeModel)
%   schedule  - initial proposed injection schedule (which also includes
%               information on the placement of wells)
%   min_wvals - minimum allowed well control values during injection phase
%   max_wvals - maximum allowed well control values during injection phase
%   varargin  - An optional number of paired arguments on the form: 
%               'option', value.   Possible options are:
%               - dryrun:       if 'true', no optimization will take place (only
%                               the associated data structures will be set up)
%               - obj_scaling:  scaling factor of the objective function (if
%                               left empty, will be computed internally)
%               - leak_penalty: penalty leak factor (default: 10)
%               - last_control_is_migration: (true or false)
%                               if 'true', the last control will be
%                               constrained to zero rate (assumed to
%                               represent the migration phase)
%               - obj_fun:      specify the objective function to maximise.  If
%                               left empty, a default objective function will
%                               be used that aims to maximise injected CO2
%                               while minimizing leakage.  The objective
%                               function, if specified, should take the
%                               following arguments: 
%                               - wellSols: vector of well solutions for all
%                                           simulation timesteps
%                               - states:   vector of simulation states for
%                                           all timesteps
%                               - schedule: the injection schedule
%                               - varargin: optional arguments should support
%                                           'ComputePartials', and 'tStep'.
%                                           If 'ComputePartials' is true,
%                                           computed values will be ADI
%                                           (automatic differentiation
%                                           objects). 
%                               The function shall return a cell array
%                               containing the objective value for each timestep.
%
% RETURNS:
%   optim   - Data structure containing all information about the final,
%             optimized scenario (including wells and the optimized schedule).
%   init    - Data structure containing all information about the starting
%             point scenario (before optimization of rates).
%   history - Information about the rate optimization process
% 
% EXAMPLE:
%   For an example, refer to the sample script 'optimizeUtsira'.
%{
Copyright 2009-2022 SINTEF Digital, Mathematics & Cybernetics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}

   
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
   opt.lineSearchMaxIt  = 10;
   opt.gradTol          = 1e-3;
   opt.objChangeTol     = 5e-4;
   opt.obj_fun = @(dummy) 0;
   opt = merge_options(opt, varargin{:});

   opt.obj_fun = @(wellSols, states, schedule, varargin) ...
                  leak_penalizer(model, wellSols, states, schedule, opt.leak_penalty, ...
                                 varargin{:});
   
   opt = merge_options(opt, varargin{:});
   
   num_wells = numel(schedule.control(1).W);
   assert(numel(max_wvals) == num_wells);
   assert(numel(min_wvals) == num_wells);
   
   init.schedule = schedule;
   
   % Check that schedule is in line with assumptions regarding well
   % control signs and types.
   assert( sameControlSigns( schedule, 1:numel(schedule.control) ), ...
       'Each well control SIGN must remain constant over all control periods.' )
   
   if opt.last_control_is_migration
       assert( all(strcmpi({schedule.control(end).W.type},'rate')), ...
            'All well control TYPES must be "rate" during the migration period.' )
       inj_control_periods = 1:numel(schedule.control)-1;
   else
       inj_control_periods = 1:numel(schedule.control);
   end
   assert( sameControlTypes( schedule, inj_control_periods ), ...
            'A well control TYPE cannot change during the injection phase.' )

   % Compute initial objective value (if required for scaling)
   if isempty(opt.obj_scaling)
      [init.wellSols, init.states] = ...
          simulateScheduleAD(initState, model, schedule, ...
                             'NonLinearSolver', NonLinearSolver('useRelaxation', true));
      
      init.obj_val_steps = cell2mat(opt.obj_fun(init.wellSols, ...
                                                init.states, ...
                                                init.schedule));
      init.obj_val_total = sum(init.obj_val_steps);
      
      % Use value of objective function before optimization as scaling.
      opt.obj_scaling = abs(init.obj_val_total);
   end
   
   % Return if optimization should be skipped
   if opt.dryrun
      optim = init;
      history = [];
      return;
   end
   
   % Define limits, scaling and objective function

   scaling.boxLims = [min_wvals(:), max_wvals(:)];
   scaling.obj     = opt.obj_scaling;
   scaling_mig     = [];
   
   if opt.last_control_is_migration
       % We ensure the min box limit applied to the well controls during
       % this last control period is "zero". As such, any wells that were
       % 'bhp' controlled during the injection period(s) are assigned a
       % "zero" min rate. (The max box limit could in theory be any value,
       % since rates will be constrained to their min values anyways.) Any
       % wells that are producers are assigned with negative rate limits.
       scaling_mig.boxLims = scaling.boxLims;
       scaling_mig.boxLims(:,1) = sqrt(eps); % @@ required?
       bhpTypeInx = strcmpi({schedule.control(1).W.type},'bhp');
       scaling_mig.boxLims( bhpTypeInx, 1 ) = sqrt(eps);
       scaling_mig.boxLims( bhpTypeInx, 2 ) = 2*sqrt(eps);
       scaling_mig.boxLims = [schedule.control(1).W.sign]' .* scaling_mig.boxLims;
   end
   
   obj_evaluator = @(u) evaluate_objective(u, opt.obj_fun, model, initState, ...
                                           schedule, scaling, scaling_mig, ...
                                           opt.last_control_is_migration ); 

   % Define constraints
   
   linEqS = [];
   if opt.last_control_is_migration
      % Constrain rates of last step to "zero" (i.e., their minimum values)
      min_rates = scaling_mig.boxLims(:,1);
      linEq = struct('A', eye(num_wells), 'b', min_rates(:));
      linEqS = setupConstraints(linEq, schedule, scaling_mig);

      % keep only the relations pertaining to the last control step
      last_step_ix = num_wells * (numel(schedule.control) - 1) + 1;
      linEqS.A = linEqS.A(last_step_ix:end, :);
      linEqS.b = linEqS.b(last_step_ix:end);
   end
   
   
   % Call optimization routine
   
   u = schedule2control(schedule, scaling);
   if opt.last_control_is_migration
       % We scale the well control values that correspond to the migration
       % period using scaling_mig.
       sch_mig.control = schedule.control(end);
       u(end-num_wells+1:end) = schedule2control(sch_mig, scaling_mig);
   end
   [~, u_opt, history] = unitBoxBFGS(u, obj_evaluator, 'linEq', linEqS, ...
                                     'lineSearchMaxIt', opt.lineSearchMaxIt, ...
                                     'gradTol',         opt.gradTol, ...
                                     'objChangeTol',    opt.objChangeTol);
                                     %'lineSearchMaxIt', 10, 'gradTol', 2e-3);
   
   % Preparing solution structures
   
   optim.schedule = control2schedule(u_opt, schedule, scaling);
   if opt.last_control_is_migration
       % We scale the well control values that correspond to the migration
       % period using scaling_mig.
       sch_mig.control = schedule.control(end);
       sch_mig = control2schedule(u_opt(end-num_wells+1:end), sch_mig, scaling_mig);
       optim.schedule.control(end) = sch_mig.control;
   end
   
   [optim.wellSols, optim.states] = ...
       simulateScheduleAD(initState, ...
                          model, ...
                          optim.schedule, ...
                          'NonLinearSolver', NonLinearSolver('useRelaxation', true));
   
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
   tot_inj = 0; % @@
   for step = 1:numSteps
      sol = wellSols{tSteps(step)};
      state = states{tSteps(step)}; %@@ +1?
      nW = numel(sol);
      pBHP = zeros(nW, 1); % place holder
      qGs = vertcat(sol.qGs);
      qWs = vertcat(sol.qWs);
      p = state.pressure;
      sG = state.s(:,2);
      sGmax = state.sGmax;
      if opt.ComputePartials
         %[p, sG, qWs, qGs, pBHP] = initVariablesADI(p, sG, qWs, qGs, pBHP);%#ok
         [p, sG, sGmax, qWs, qGs, pBHP] = initVariablesADI(p, sG, sGmax, qWs, qGs, pBHP);%#ok
      end
      dt = dts(step);
      injInx = (vertcat(sol.sign) > 0);
      obj{step} = dt * spones(ones(1, nW)) * ((1-penalty) * injInx .* qGs);

      tot_inj = tot_inj +  dt * spones(ones(1, nW)) * ( injInx .* qGs);
      
      if (tSteps(step) == num_timesteps)
         bG = model.fluid.bG(p);
         if ~isfield(model.rock,'ntg')
             model.rock.ntg = ones(model.G.cells.num,1); % in case ntg doesn't exist
         end
         pvol = model.G.cells.volumes .* model.G.cells.H .* model.rock.poro .* model.rock.ntg;
         vol = ones(1, model.G.cells.num) * (pvol .* model.fluid.pvMultR(p) .* bG .* sG);
         obj{step} = obj{step} + penalty * vol;
         if ~opt.ComputePartials
            fprintf('Total injected: %f\n', double(tot_inj));
            fprintf('Total leaked: %f\n', double(tot_inj - vol));
            fprintf('Score: %f\n\n', double(tot_inj) - penalty * double(tot_inj-vol));
         end
      end
      obj{step} = obj{step} * model.fluid.rhoGS/1e12; % Gt CO2
   end
end

% ----------------------------------------------------------------------------

function [val, der, wellSols, states] = ...
       evaluate_objective(u, obj_fun, model, initState, schedule, scaling, ...
                          scaling_mig, last_control_is_migration) 
   
% With possiblity of separate scalings applied to the injection period(s)
% and the migration period. This is required when a well is bhp-controlled
% during the injection period(s), and then rate-controlled during the
% migration period.

   if last_control_is_migration
       assert( ~isempty(scaling_mig), ...
           'Need to supply separate scaling for migration period.' )
       num_wells            = numel(schedule.control(1).W);
       schedule_mig.control = schedule.control(end);
   end

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
   if last_control_is_migration
       % We scale the well control values that correspond to the migration
       % period using scaling_mig.
       sch_mig = control2schedule(u(end-num_wells+1:end), schedule_mig, scaling_mig);
       schedule.control(end) = sch_mig.control;
       clear sch_mig
   end
   
   % run simulation:
   [wellSols, states] = ...
       simulateScheduleAD(initState, model, schedule, ...
                          'NonLinearSolver', NonLinearSolver('useRelaxation', true));
   
   % compute objective:
   vals = obj_fun(wellSols, states, schedule);
   val  = sum(cell2mat(vals))/abs(objScaling);
   val  = full(val); % ensure 'val' is a regular scalar, not a sparse one
   
   % run adjoint:
   if nargout > 1
      objh = @(tstep,varargin)obj_fun(wellSols, states, schedule, 'ComputePartials', true, 'tStep', tstep);
      g    = computeGradientAdjointAD(initState, states, model, schedule, objh);
      % scale gradient:
      der = scaleGradient(g, schedule, boxLims, objScaling);
      der = vertcat(der{:});
      if last_control_is_migration
          schedule_mig.control  = schedule.control(end);
          boxLims_mig           = scaling_mig.boxLims;
          der_mig     = scaleGradient({g{end}}, schedule_mig, boxLims_mig, objScaling);
          der         = vertcat(der(1:end-num_wells), der_mig{:});
      end
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

function flag = sameControlTypes( schedule, control_periods )

    % Output is true/false depending on whether or not an individual well
    % is of the same control type (i.e., 'bhp' or 'rate') over the
    % specified control periods.
    
    % Example: control_period = [1 2], or [1 2 3], or [2 3], etc.
    
    assert( max(control_periods) <= numel(schedule.control), ...
        'The schedule does not contain control period %d.', max(control_periods) ) 
    
    nc = numel(control_periods);
    tmp = [];
    for i=1:nc-1
        a = { schedule.control( control_periods(i) ).W.type };
        b = { schedule.control( control_periods(i+1) ).W.type };
        tmp = [tmp, isequal(a, b)];
    end
    flag = all(tmp); % for flag to be true, all values of tmp must be 1 (true).

end

% ----------------------------------------------------------------------------

function flag = sameControlSigns( schedule, control_periods )

    % Output is true/false depending on whether or not an individual well
    % is of the same control sign (i.e., +1 or -1) over the specified
    % control periods.
    
    % Example: control_period = [1 2], or [1 2 3], or [2 3], etc.

    assert( max(control_periods) <= numel(schedule.control), ...
        'The schedule does not contain control period %d.', max(control_periods) ) 

    nc = numel(control_periods);
    tmp = [];
    for i=1:nc-1
        a = { schedule.control( control_periods(i) ).W.sign };
        b = { schedule.control( control_periods(i+1) ).W.sign };
        tmp = [tmp, isequal(a, b)];
    end
    flag = all(tmp); % for flag to be true, all values of tmp must be 1 (true).

end

