function obj = pressurePenalizer(model, states, schedule, penalty, plim, varargin)
% Penalize pressure of specific or all cells.
%
% SYNOPSIS:
%   obj = pressurePenalizer(model, states, schedule, penalty, plim)
%   obj = pressurePenalizer(model, states, schedule, penalty, plim, 'cells', cells)
%
% DESCRIPTION:
%   Evaluate an objective function which penalizes specific or all cell
%   pressures that surpasses their pressure limit, plim. Can also penalize
%   bottom-hole pressure of each well if it surpasses their pressure limit,
%   BHPlim.
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
%   penalizeBHP         - true or false flag used to determine whether or
%                       not to penalize the bottom-hole pressure of each
%                       well. Default is false.
%   BHPlim              - the limit for the bottom-hole pressure of each
%                       well. Required if penalizeBHP is true. Default is
%                       empty.
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

% Format of objective function:
%   obj = max(0, sign(p - plim)) * penalty * (p - plim)^k
% obj is computed for each or for a specific time step

   opt.cells            = [];
   opt.penalizeBHP      = false;
   opt.BHPlim           = [];       % only used if penalizeBHP is true
   opt.ComputePartials  = false;
   opt.tStep            = [];
   opt = merge_options(opt, varargin{:});
   
   if opt.penalizeBHP
      assert( ~isempty(opt.BHPlim), ...
          'In order to penalize the BHP of the wells, the BHP limits must be specified.' )
      assert( numel(opt.BHPlim) == numel(schedule.control(1).W), ...
          'Expected %d BHP limits for %d wells', numel(schedule.control(1).W), ...
          numel(schedule.control(1).W) ) 
      assert( all(isfinite(opt.BHPlim)), ...
          'Bottom-hole pressure limits must be finite values.' )
   end
   
   if isempty(opt.cells)
      cells = 1:model.G.cells.num;
   else
      cells = opt.cells;
      plim  = plim(cells);
   end
   
   assert( all(isfinite(plim)), ...
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
      p     = state.pressure;       % size G.cells.num x 1 (i.e., all grid cell pressures)
      if opt.penalizeBHP
        pBHP  = [state.wellSol.bhp]'; % size number of wells x 1
      end
      if opt.ComputePartials
        sG = state.s(:,2);   % place holders
        sGmax = state.sGmax; % place holders
        nW = numel(schedule.control(1).W);
        if ~opt.penalizeBHP
            pBHP = zeros(nW, 1); % place holders
        end
        qGs = zeros(nW, 1);  % place holders
        qWs = qGs;           % place holders
        
        % Depending on the model, or more specifically the equation
        % variables, we include sGmax as an ADI variable. There is likely a
        % better way to perform these conditional statements, rather than
        % checking the model class.
        if strcmpi(class(model),'CO2VEBlackOilTypeModel')
            [p, ~, ~, ~, ~, ~] = initVariablesADI(p, sG, sGmax, qWs, qGs, pBHP); % size of p is G.cells.num (i.e., all grid cells)
        elseif strcmpi(class(model),'WaterGasMultiVEModel')
            if opt.penalizeBHP
                [p, ~, ~, ~, pBHP] = initVariablesADI(p, sG, qWs, qGs, pBHP); % size of p is G.cells.num (i.e., all grid cells), size of pBHP is number of wells
            else
                [p, ~, ~, ~, ~] = initVariablesADI(p, sG, qWs, qGs, pBHP); % size of p is G.cells.num (i.e., all grid cells)
            end
        else
            error(['Did not expect the model class you are using. ',...
                'Modify conditional statements in pressurePenalizer as needed.'])
        end
      end
      A = sparse(1:numel(cells), cells, 1, numel(cells), model.G.cells.num);
      p = A * p; % p.val is now for specified cells, and the size of the
                 % jacobians (if we are computing partials) is num
                 % specified cells x number of grid cells
      
      % First penalize cell pressures
      tmp = max(0, sign(p - plim)) .* penalty .* (p - plim).^k/1e12; % scaling to MPa
      % tmp will be numeric if all values of tmp are computed to be 0 (tmp
      % will still be an adi if at least one value of tmp is computed to be
      % non-zero). Thus, we ensure tmp gets converted to an adi with its
      % expected zero derivatives.
      if opt.ComputePartials && isnumeric(tmp)
         tmp = double2ADI(tmp, p);
      end
      tmp = tmp .* model.G.cells.volumes(cells);
      obj{step} = sum( tmp )./sum(model.G.cells.volumes(cells));
      
      % Then penalize bottom-hole pressures (if applicable) and add to obj
      % function value
      if opt.penalizeBHP
          clear tmp
          tmp  = max(0, sign(pBHP - opt.BHPlim)) .* penalty .* (pBHP - opt.BHPlim).^k/1e12;
          if opt.ComputePartials && isnumeric(tmp)
              tmp = double2ADI(tmp, pBHP);
          end
          % tmp doesn't need to be normalized by cell volumes. Simply sum.
          obj{step} = sum( tmp ) + obj{step};
      end
      
      dt = dts(step);
      obj{step} = obj{step} * dt; % obj jacobians correspond to the grid size.

      
      % Keep track of amount over or amount under plim at each time step:
      % (Used for reporting purposes only)
      if opt.penalizeBHP
          P     = [p; pBHP];
          PLIM  = [plim; opt.BHPlim]; % use different variable name to avoid overwriting plim
      else
          [P, PLIM] = deal(p, plim);
      end
      PD = value(P); % we take double here to get cinx of possible ADI
      max_amount_surp(step) = max(0,max(PD-PLIM));
      if max_amount_surp(step) > 0
          [~,cinx] = max(PD-PLIM);
          max_amount_surp_cinx(step) = cinx; % this cinx does not correspond
                                             % to the grid's cell index,
                                             % rather it is the index of p
                                             % that surpassed its plim the
                                             % most
      end
      min_amount_under(step) = max(0,min(PLIM-PD));
      if min_amount_under(step) > 0
         [~,cinx] = min(max(0,(PLIM-PD)));
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
                 msg1 = 100*val/PLIM(cinx);
             else
                 [val,tinx] = min(min_amount_under);
                 cinx = min_amount_under_cinx(tinx);
                 msg2 = 100*val/PLIM(cinx);
             end
             fprintf('Surpassed Plimit by %f (percent) of Plimit.\n', msg1)
             fprintf('Approached Plimit by %f (percent) of Plimit.\n', msg2)
         end
      end

   end
end
