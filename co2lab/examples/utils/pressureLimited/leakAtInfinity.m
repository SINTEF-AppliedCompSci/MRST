function obj = leakAtInfinity(model, wellSols, states, schedule, ...
    penalty, surf_press, rho_water, ta, varargin)
% computes mass (in Gt) leaked by time infinity as long as penalty is set
% to 1.
% NB: obj value is in terms of a negative mass...
%
% pressure is assumed to be hydrostatic at time infinity.

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

   opt.ComputePartials = false;
   opt.tStep = [];
   opt = merge_options(opt, varargin{:});
   
   assert(penalty == 1, 'Expected penalty value of 1.')
   
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
   
   if ~isfield(model.rock,'ntg')
      model.rock.ntg = ones(model.G.cells.num,1); % in case ntg doesn't exist
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
      p_future = rho_water * norm(gravity()) * model.G.cells.z + surf_press;
      sG = state.s(:,2);
      sF = state.s(:,1);
      sGmax = state.sGmax;
      rs = zeros(model.G.cells.num, 1); % place holder
      if isfield(model.fluid,'rsSat') % i.e., dissolution is true
        rs = state.rs;
        warning('MassAtInfinity not well tested for solubility trapping.')
      end
      if opt.ComputePartials
         [p, sG, sGmax, qWs, qGs, pBHP] = initVariablesADI(p, sG, sGmax, qWs, qGs, pBHP);%#ok
      end
      dt = dts(step);
      injInx = (vertcat(sol.sign) > 0);
      
      % calculate the "(1 - C) M^inj" part of obj fun J
      obj{step} = dt * spones(ones(1, nW)) * ((-penalty) * injInx .* qGs);

      krull = krull +  dt * spones(ones(1, nW)) * ( injInx .* qGs);
      
      if (tSteps(step) == num_timesteps)
          
         % calculate the "C M^a" part of obj fun J, where M^a is the mass
         % accumulated (or remaining) at t = infinity. To determine this
         % mass, we calculate the volume expected to be retained using
         % spill-pt dynamics.
         dh = [];   % does not consider subscale trapping
         will_stay = massAtInfinity( model.G, model.rock, p, sG, sGmax, sF, rs, ...
             model.fluid, ta, dh, 'p_future', p_future ); % kg
         vol_inf   = will_stay / model.fluid.rhoGS;       % m3, ref. depth
         
         % J = (1-C) M^i + C M^a
         obj{step} = obj{step} +  penalty*vol_inf;
         if ~opt.ComputePartials
            fprintf('Total injected: %f (m3)\n', double(krull));
            %fprintf('Total leaked (by infinity): %f (m3)\n', double(krull - vol_inf));
            %fprintf('Penalty: %f \n', penalty)
            fprintf('Score: %f (m3)\n\n', double(krull) - penalty * double(krull-vol_inf));
         end
      end
      obj{step} = obj{step} * model.fluid.rhoGS/1e12; % vol * rho, in Gt.
   end
end
