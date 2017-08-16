function [optim, init, history] = optimizeRates(initState, model, schedule, ...
                                                min_wvals, max_wvals, varargin)
%
% Compute an optimal set of well controls ('rate' or 'bhp') for a proposed
% injection/migration scenario.
% 
% SYNOPSIS:
%   function [optim, init, history] = ...
%      optimizeRates(initState, model, schedule, min_wvals, max_wvals, varargin)
%
% DESCRIPTION:
%
% PARAMETERS:
%   initState - initial state
%   model     - simulation model (instance of CO2VEBlackOilTypeModel)
%   schedule  - initial proposed injection schedule (which also includes
%               information on the placement of wells)
%   min_wvals - minimum allowed well control values during injection period
%   max_wvals - maximum allowed well control values during injection period
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
Copyright 2009-2017 SINTEF ICT, Applied Mathematics.

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
   % NB: if well control types for the injection and migration periods are
   % the same, the same boxLims will be applied to the well controls during
   % both control periods. However, if the well control types are different
   % (e.g., 'bhp' during injection/production, then 'rate' during
   % migration), we require two separate scalings. We assume the well
   % controls during the migration period are type 'rate'.

   scaling.boxLims = [min_wvals(:), max_wvals(:)];
   scaling.obj     = opt.obj_scaling;
   
   assert( numel(schedule.control)==2 , 'Expected two control periods.' )
   assert( all(strcmpi({schedule.control(2).W.type},'rate')), ...
       'Expected well controls during migration period to be type "rate".' )
   
   % Flag for whether or not control types are the same in both periods.
   sameControlTypes = all(strcmpi({schedule.control(1).W.type}, ...
       {schedule.control(2).W.type}));

   if sameControlTypes
       scaling_mig = [];
   else
       scaling_mig.boxLims = scaling.boxLims;
       bhpTypeInx = strcmpi({schedule.control(1).W.type},'bhp');
       scaling_mig.boxLims( bhpTypeInx, 1 ) = -sqrt(eps);
       scaling_mig.boxLims( bhpTypeInx, 2 ) = -2*sqrt(eps);
   end
   
   obj_evaluator = @(u) evaluate_objective(u, opt.obj_fun, model, initState, ...
                                           schedule, scaling, scaling_mig); 

   % Define constraints
   
   linEqS = [];
   if opt.last_control_is_migration
      % Constrain rates of last step to zero
      if sameControlTypes
          min_rates = scaling.boxLims(:,1);
      else
          min_rates = scaling_mig.boxLims(:,1);
      end
      linEq = struct('A', eye(num_wells), 'b', min_rates(:));
      if sameControlTypes
          linEqS = setupConstraints(linEq, schedule, scaling);
      else
          linEqS = setupConstraints(linEq, schedule, scaling_mig);
      end
      
      
      % keep only the relations pertaining to the last control step
      last_step_ix = num_wells * (numel(schedule.control) - 1) + 1;
      linEqS.A = linEqS.A(last_step_ix:end, :);
      linEqS.b = linEqS.b(last_step_ix:end);
   end
   
   
   % Call optimization routine
   
   if sameControlTypes
       u = schedule2control(schedule, scaling);
   else
       sch1.control = schedule.control(1);
       sch2.control = schedule.control(2);
       u_1 = schedule2control(sch1, scaling);
       u_2 = schedule2control(sch2, scaling_mig);
       u = [u_1; u_2];
   end
   [~, u_opt, history] = unitBoxBFGS(u, obj_evaluator, 'linEq', linEqS, ...
                                     'lineSearchMaxIt', opt.lineSearchMaxIt, ...
                                     'gradTol',         opt.gradTol, ...
                                     'objChangeTol',    opt.objChangeTol);
                                     %'lineSearchMaxIt', 10, 'gradTol', 2e-3);
   
   % Preparing solution structures
   
   if sameControlTypes
       optim.schedule = control2schedule(u_opt, schedule, scaling);
   else
       sch1 = control2schedule(u_opt(1:end/2), sch1, scaling);
       sch2 = control2schedule(u_opt(end/2+1:end), sch2, scaling_mig);
       optim.schedule.control(1) = sch1.control;
       optim.schedule.control(2) = sch2.control;
       optim.schedule.step = schedule.step;
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
      obj{step} = obj{step} * model.fluid.rhoGS/1e12; %@@
   end
end

% ----------------------------------------------------------------------------

function [val, der, wellSols, states] = ...
       evaluate_objective(u, obj_fun, model, initState, schedule, scaling, scaling_mig) 
   
% With possiblity of separate scalings applied to two different control
% periods. This is required when a well is bhp-controlled during the
% injection period, and is rate-controlled during the migration period.

   if numel(schedule.control)==1
       % Since there is only one control period, there is no difference
       % between control types.
       sameControlTypes = true;
   
   elseif numel(schedule.control)==2
       % We check whether there is a difference between well control types
       % between the first (injection) and second (migration) control
       % periods.
       sameControlTypes = all(strcmpi({schedule.control(1).W.type}, {schedule.control(2).W.type}));
       assert( ~isempty(scaling_mig), 'Need non-empty separate scaling for migration period.' )
       
   else
       error( 'Not yet implemented for more than 2 control periods.' )
   end
   
   if ~sameControlTypes
       % separate the two controls
       sch1.control = schedule.control(1);
       sch2.control = schedule.control(2);
       % where u = [u1; u2]
   end

   minu = min(u);
   maxu = max(u);
   if or(minu < -eps , maxu > 1+eps)
      warning('Controls are expected to lie in [0 1]\n')
   end
   
   if sameControlTypes
       boxLims = scaling.boxLims;
   else
       boxLims1 = scaling.boxLims;
       boxLims2 = scaling_mig.boxLims;
   end

   if isfield(scaling, 'obj')
      objScaling = abs(scaling.obj);
   else
      objScaling = 1;
   end
   
   % update schedule:
   if sameControlTypes
       schedule = control2schedule(u, schedule, scaling);
   else
       % the first half of the well controls u belongs to the first control
       % period, and the second half belongs to the second control period.
       % We scale these well controls separately.
       sch1 = control2schedule(u(1:end/2), sch1, scaling);
       sch2 = control2schedule(u(end/2+1:end), sch2, scaling_mig);
       schedule.control(1) = sch1.control;
       schedule.control(2) = sch2.control;
   end
   
   % run simulation:
   [wellSols, states] = ...
       simulateScheduleAD(initState, model, schedule, ...
                          'NonLinearSolver', NonLinearSolver('useRelaxation', true));
   
   % compute objective:
   vals = obj_fun(wellSols, states, schedule);
   val  = sum(cell2mat(vals))/abs(objScaling);

   % run adjoint:
   if nargout > 1
      objh = @(tstep)obj_fun(wellSols, states, schedule, 'ComputePartials', true, 'tStep', tstep);
      g    = computeGradientAdjointAD(initState, states, model, schedule, objh);
      % scale gradient:
      if sameControlTypes
          der = scaleGradient(g, schedule, boxLims, objScaling);
          der = vertcat(der{:});
      else
          der1 = scaleGradient({g{1}}, sch1, boxLims1, objScaling);
          der2 = scaleGradient({g{2}}, sch2, boxLims2, objScaling);
          der = vertcat(der1{:}, der2{:});
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
    

