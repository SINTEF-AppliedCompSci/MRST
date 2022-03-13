function [obj, Mi_tot, Ma] = leak_penalizer_at_infinity_Rerun(model, wellSols, states, schedule, ...
    penalty, surf_press, rho_water, ta, varargin)
% Purpose: to compute the obj value at each time-step:
% J(t) = Mi * (1-C) + C Ma_inf
%
% Output 'obj' not designed to be used to evaluate objective function
% (i.e., taking sum(obj) is not computing the total objective function).
%
% Pressure is assumed to be hydrostatic at time infinity.

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

   opt.tStep = [];
   opt.plotsOn = false;
   opt.h = 48; % figure handle for plotting
   opt.dh = [];
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
   [Mi_tot, Ma] = deal(zeros(numSteps,1));
   krull = 0; % @@
   
   if ~isfield(model.rock,'ntg')
      model.rock.ntg = ones(model.G.cells.num,1); % in case ntg doesn't exist
   end 
   
   for step = 1:numSteps
      sol = wellSols{tSteps(step)};
      state = states{tSteps(step)}; %@@ +1?
      nW = numel(sol);
      qGs = vertcat(sol.qGs);
      p = state.pressure;
      p_future = rho_water * norm(gravity()) * model.G.cells.z + surf_press;
      sG = state.s(:,2);
      sF = state.s(:,1);
      sGmax = state.sGmax;
      rs = zeros(model.G.cells.num, 1); % place holder
      if isfield(model.fluid,'rsSat') % i.e., dissolution is true
         rs = state.rs;
      end
      
      % J = Mi * (1-C):
      dt = dts(step);
      injInx = (vertcat(sol.sign) > 0);
      obj{step} = dt * spones(ones(1, nW)) * ((1-penalty) * injInx .* qGs);
      
      % Mass injected so far:
      Mi_tot(step) = dt * spones(ones(1, nW)) * (injInx .* qGs) * model.fluid.rhoGS/1e12; % Gt
      if step>1
        Mi_tot(step) = Mi_tot(step - 1) + Mi_tot(step); % total injected, Gt
      end

      krull = krull +  dt * spones(ones(1, nW)) * ( injInx .* qGs);
      
      % J = J + C * Ma_inf:
      % where M^a_inf is the mass accumulated (or remaining) at t = infinity.
      will_stay = massAtInfinity( model.G, model.rock, p, sG, sGmax, sF, rs, ...
          model.fluid, ta, opt.dh, 'p_future', p_future );  % kg
      vol_inf   = will_stay / model.fluid.rhoGS;            % m3, ref. depth
      
      obj{step} = obj{step} + penalty * vol_inf;
      Ma(step) = vol_inf * model.fluid.rhoGS/1e12; % Gt
      
      if (tSteps(step) == num_timesteps)
         fprintf('Total injected: %f (m3)\n', double(krull));
         fprintf('Total leaked (by infinity): %f (m3)\n', double(krull - vol_inf));
         fprintf('Penalty: %f \n', penalty)
         fprintf('Score: %f (m3)\n\n', double(krull) - penalty * double(krull-vol_inf));
         % to make final plot
         [~] = massAtInfinity( model.G, model.rock, p, sG, sGmax, sF, rs, ...
          model.fluid, ta, opt.dh, 'p_future', p_future, 'plotsOn', opt.plotsOn, 'h',opt.h );  % kg
      end
      obj{step} = obj{step} * model.fluid.rhoGS/1e12; % vol * rho, in Gt.
      
   end

end
