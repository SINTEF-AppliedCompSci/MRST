function obj = maxPressureViolation(model, states, schedule, plim, varargin)
% states.pressure is a cell array.
% schedule is only used for time steps.
% plim is a cell array.

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

% format of objective function:
%   obj = [... 0 0 0 dp 0 0 0 0 0 ...]
% where dp is max pressure violation (in Pascals)

% obj is computed for each time step (for both inj and mig periods), but is
% only non-zero at time step where max pressure violation occurred.

% dp can be negative (indicating plim was not surpassed)

   
   opt.ComputePartials  = false;
   opt.tStep            = [];
   opt.cells            = [1:model.G.cells.num]';
   opt.hfig             = [];  % figure handle (a scalar numeric value) for
                               % plotting location of max pressure
                               % violation. If empty, no plot generated.
                               % Default is empty.
   opt.tsteps           = ones(numel(states),1);
   
   opt = merge_options(opt, varargin{:});
   
   num_timesteps = numel(schedule.step.val);
   
   dpmaxstep=-inf(num_timesteps,1);
   for i=1:num_timesteps
      if(opt.tsteps(i)==1) 
        dpmaxstep(i)=max(states{i}.pressure(opt.cells)-plim(opt.cells)); 
      end
   end
   [dpmaxa,dmaxstep]=max(dpmaxstep);
   
   tSteps = opt.tStep;
   if isempty(tSteps)
      numSteps = numel(states);
      tSteps = (1:numSteps)';
      %dts = schedule.step.val;
   else
      assert(numel(tSteps) == 1);
      numSteps = 1;
      %dts = schedule.step.val(opt.tStep);
   end
   
   obj = repmat({[]}, numSteps, 1);
  
   max_amount_surp      = zeros(numSteps, 1);
   max_amount_surp_cinx = zeros(numSteps, 1);
   min_amount_under     = zeros(numSteps, 1);
   min_amount_under_cinx= zeros(numSteps, 1);
   for step = 1:numSteps     
      state = states{tSteps(step)}; %@@ +1?      
      p = state.pressure;
      % keep track of amount over or amount under plim at each time step
      max_amount_surp(step) = max(0,max(p(opt.cells)-plim(opt.cells)));
      if max_amount_surp(step) > 0
          [~,cinx] = max(p(opt.cells)-plim(opt.cells));
          max_amount_surp_cinx(step) = opt.cells(cinx);
      end
      min_amount_under(step) = max(0,min(plim(opt.cells)-p(opt.cells)));
      if min_amount_under(step) > 0
         [~,cinx] = min(max(0,(plim(opt.cells)-p(opt.cells))));
         min_amount_under_cinx(step) = opt.cells(cinx);
      end
      if opt.ComputePartials
        sG = state.s(:,2);   % place holders
        sGmax = state.sGmax; % place holders
        nW = numel(schedule.control(1).W);
        pBHP = zeros(nW, 1); % place holders
        qGs = pBHP;          % place holders
        qWs = pBHP;          % place holders
        [p, ~, ~, ~, ~, ~] = initVariablesADI(p, sG, sGmax, qWs, qGs, pBHP); 
      end
      
      tmp=max(p(opt.cells)-plim(opt.cells)); % Pascals
      if(tSteps(step)==dmaxstep)
        assert(double(tmp)==dpmaxa);
        obj{step} = tmp;
      else
        obj{step}= 0*tmp;  
      end
      
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
             fprintf('   in cell %d    at timestep %d \n', cinx, dmaxstep)
             
             % plot grid location where max pressure violation occurred
             % (optional)
             if ~isempty(opt.hfig)
                 if ~ishandle(opt.hfig)
                    figure(opt.hfig)
                 else
                    % Avoid stealing focus if figure already exists
                    set(0, 'CurrentFigure', opt.hfig);
                    clf
                 end
                 %figure(min(opt.cells)), clf
                 plotGrid(model.G, 'facecolor','none', 'edgealpha',0.1),
                 plotCellData(model.G, model.G.cells.z, opt.cells, ...
                     'facecolor','y', 'edgealpha',0.1) % the cells being checked
                 if any(max_amount_surp > 0)
                    plotCellData(model.G, model.G.cells.z, cinx, 'facecolor','r')
                    title(['cell ',num2str(cinx),', timestep ',num2str(dmaxstep)])
                 else
                    title('pressure limit not violated in highlighted region')
                 end
             end
             
         end
      end

   end
end
