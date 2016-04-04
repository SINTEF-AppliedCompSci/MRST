function [obj, Mi_tot, Ma] = leak_penalizer_at_infinity_Rerun(model, wellSols, states, schedule, ...
    penalty, surf_press, rho_water, varargin)
% copy of local helper function from optimizeFormation:

% Purpose: to compute the obj value at each time-step:
% J(t) = Mi * (1-C) + C Ma_inf

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
   [Mi_tot, Ma] = deal(zeros(numSteps,1));
   krull = 0; % @@
   
   % compute ta here, then use to pass into vol_at_infinity so it is not
   % re-computed each time vol_at_infinity is called.
   ta = trapAnalysis(model.G, false);
   
   for step = 1:numSteps
      sol = wellSols{tSteps(step)};
      state = states{tSteps(step)}; %@@ +1?
      nW = numel(sol);
      pBHP = zeros(nW, 1); % place holder
      qGs = vertcat(sol.qGs);
      qWs = vertcat(sol.qWs);
      %p = state.pressure;
      p = compute_hydrostatic_pressure(model.G, rho_water, surf_press);
      sG = state.s(:,2);
      if opt.ComputePartials
         %[p, sG, pBHP, qWs, qGs] = initVariablesADI(p, sG, pBHP, qWs, qGs);%#ok
         [p, sG, qWs, qGs, pBHP] = initVariablesADI(p, sG, qWs, qGs, pBHP);%#ok
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
      bG = model.fluid.bG(p);
      if ~isfield(model.rock,'ntg')
             model.rock.ntg = ones(model.G.cells.num,1); % in case ntg doesn't exist
         end 
      %pvol = model.G.cells.volumes .* model.G.cells.H .* model.rock.poro .* model.rock.ntg;
      %vol = ones(1, model.G.cells.num) * (pvol .* model.fluid.pvMultR(p) .* bG .* sG); % @@ what pressure?
      vol_inf = vol_at_infinity( model.G, ...
             model.rock, ... % or model.rock.poro .* model.fluid.pvMultR(p), ...
             sG .* bG, ...
             model.fluid.res_water, model.fluid.res_gas, 'ta',ta); % using possible ADI variables
      obj{step} = obj{step} + penalty * vol_inf;
      Ma(step) = vol_inf * model.fluid.rhoGS/1e12; % Gt
      
      if (tSteps(step) == num_timesteps)
         if ~opt.ComputePartials
            fprintf('Total injected: %f (m3)\n', double(krull));
            fprintf('Total leaked: %f (m3)\n', double(krull - vol_inf));
            fprintf('Score: %f (m3)\n\n', double(krull) - penalty * double(krull-vol_inf));
         end
      end
      obj{step} = obj{step} * model.fluid.rhoGS/1e12; % vol * rho, in Gt.
      %vol_steps(step) = vol;
      %vol_inf_steps(step) = vol_inf;
   end

end

function p = compute_hydrostatic_pressure(Gt, rho_water, surface_pressure)

    p = rho_water * norm(gravity()) * Gt.cells.z + surface_pressure;
end
    
