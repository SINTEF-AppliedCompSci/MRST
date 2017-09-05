function obj = pressurePenalizer(model, states, schedule, penalty, plim, varargin)
% Penalize pressure of specific or all cells.
%
% SYNOPSIS:
%   obj = pressurePenalizer(model, states, schedule, penalty, plim)
%   obj = pressurePenalizer(model, states, schedule, penalty, plim, 'cells', cells)
%
% DESCRIPTION:
%   Evaluate an objective function which penalizes specific or all cell
%   pressures that surpasses their pressure limit, plim.
%
% REQUIRED PARAMETERS:
%   model       - The model.
%   states      - Structure containing pressure field.
%   schedule    - Structure containing time step and well details.
%   penalty     - The pressure penalty factor. A scalar value.
%   plim        - Pressure limit for all cells (size must correspond to
%                 size of state.pressure)
%
% OPTIONAL:
%   cells               - Cell indexes of cells to be assessed in objective
%                       function. If empty, all cells are considered.
%                       Default is empty.
%   ComputePartials     - true or false flag used to determine whether or
%                       not to compute the partial derivatives of obj.
%                       function w.r.t state variables, etc. Default is
%                       false.
%   tStep               - time step at which to evaluate obj. If empty, all
%                       timesteps will be used to computed obj. Default is
%                       empty.
%
% RETURNS:
%   obj         - cell array of size Nx1 (where N is number of time steps
%               that function gets evaluated at)
%
% SEE ALSO:
%   

% Format of objective function:
%   obj = max(0, sign(p - plim)) * penalty * (p - plim)^k
% obj is computed for each or for a specific time step


   opt.cells            = [];
   opt.ComputePartials  = false;
   opt.tStep            = [];
   opt = merge_options(opt, varargin{:});
   
   if isempty(opt.cells)
      cells = 1:model.G.cells.num;
   else
      cells = opt.cells;
      plim  = plim(cells);
   end
   
   assert( all(isfinite(plim(cells))), ...
       'Pressure limit of specific or all cells must be finite values.' )
   
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

      % Evaluate objective function at current time step:
      state = states{tSteps(step)}; %@@ +1?      
      p = state.pressure; % all cells in grid
      if opt.ComputePartials
        sG = state.s(:,2);   % place holders
        sGmax = state.sGmax; % place holders
        nW = numel(schedule.control(1).W);
        pBHP = zeros(nW, 1); % place holders
        qGs = pBHP;          % place holders
        qWs = pBHP;          % place holders
        [p, ~, ~, ~, ~, ~] = initVariablesADI(p, sG, sGmax, qWs, qGs, pBHP); % this is for all cells in grid
      end
      A = sparse(1:numel(cells), cells, 1, numel(cells), model.G.cells.num);
      p = A * p; % p.val is now for specified cells, and the size of the
                 % jacobians (if we are computing partials) is num
                 % specified cells x number of grid cells
      dt = dts(step);
      tmp = max(0, sign(p - plim)) .* penalty .* (p - plim).^k/1e12; % scaling to MPa
      % tmp will be numeric if all values of tmp are computed to be 0 (tmp
      % will still be an adi if at least one value of tmp is computed to be
      % non-zero). Thus, we ensure tmp gets converted to an adi with zero
      % its expected zero derivatives.
      if opt.ComputePartials && isnumeric(tmp)
         tmp = double2ADI(tmp, p);
      end
      tmp = tmp .* model.G.cells.volumes(cells);
      obj{step} = sum( tmp )./sum(model.G.cells.volumes(cells));
      obj{step} = obj{step} * dt; % obj jacobians correspond to the grid size.

      
      % Keep track of amount over or amount under plim at each time step:
      pd = double(p); % we take double here to get cinx of possible ADI
      max_amount_surp(step) = max(0,max(pd-plim));
      if max_amount_surp(step) > 0
          [~,cinx] = max(pd-plim);
          max_amount_surp_cinx(step) = cinx; % this cinx does not correspond
                                             % to the grid's cell index,
                                             % rather it is the index of p
                                             % that surpassed its plim the
                                             % most
      end
      min_amount_under(step) = max(0,min(plim-pd));
      if min_amount_under(step) > 0
         [~,cinx] = min(max(0,(plim-pd)));
         min_amount_under_cinx(step) = cinx; % same as note above
      end
      
      
      % Reporting:
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