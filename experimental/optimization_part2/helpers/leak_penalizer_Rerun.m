function [obj, Mi, Mi_tot, Ma] = leak_penalizer_Rerun(model, wellSols, states, schedule, penalty, varargin)
% copy of local helper function from optimizeFormation:

% Purpose: to compute the obj value at each time-step:
% J(t) = Mi(t) - C (Mi(t) - Ma(t))
% J(t) = (1-C) Mi(t) + C Ma(t)
% where Mi(t) is amount injected at time t (not total injected by time t)!
% where Ma(t) is total amount in domain at time t

% Also: to return Mi and Ma at each time-step, and Mi_tot over time

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
   [Mi, Mi_tot, Ma] = deal(zeros(numSteps,1));
   krull = 0; % @@
   for step = 1:numSteps
       
      % compute full obj value at each time step
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
      
      % J = Mi(t):
      dt = dts(step);
      injInx = (vertcat(sol.sign) > 0);
      Mi(step) = dt * spones(ones(1, nW)) * (injInx .* qGs); % kg
      obj{step} = dt * spones(ones(1, nW)) * ((1-penalty) * injInx .* qGs);
      Mi_tot(step) = dt * spones(ones(1, nW)) * (injInx .* qGs) * model.fluid.rhoGS/1e12; % Gt
      if step>1
        Mi_tot(step) = Mi_tot(step - 1) + Mi_tot(step); % total injected, Gt
      end

      krull = krull +  dt * spones(ones(1, nW)) * ( injInx .* qGs);
      
      % J = J + C * Ma:
      bG = model.fluid.bG(p);
      if ~isfield(model.rock,'ntg')
            model.rock.ntg = ones(model.G.cells.num,1); % in case ntg doesn't exist
      end 
      pvol = model.G.cells.volumes .* model.G.cells.H .* model.rock.poro .* model.rock.ntg;
      vol = ones(1, model.G.cells.num) * (pvol .* model.fluid.pvMultR(p) .* bG .* sG);
      obj{step} = obj{step} + penalty * vol;
      Ma(step) = vol * model.fluid.rhoGS/1e12; % Gt
      
      if (tSteps(step) == num_timesteps)
%          bG = model.fluid.bG(p);
%          pvol = model.G.cells.volumes .* model.G.cells.H .* model.rock.poro;      
%          vol = ones(1, model.G.cells.num) * (pvol .* model.fluid.pvMultR(p) .* bG .* sG); % @@ need to account for NTG!
%          obj{step} = obj{step} + penalty * vol;
         if ~opt.ComputePartials
            fprintf('Total injected: %f (m3)\n', double(krull));
            fprintf('Total leaked: %f (m3)\n', double(krull - vol));
            fprintf('Score: %f (m3)\n\n', double(krull) - penalty * double(krull-vol));
         end
      end
      obj{step} = obj{step} * model.fluid.rhoGS/1e12; %@@ converted to Gt
      Mi(step) = Mi(step) * model.fluid.rhoGS/1e12; % Gt
   end
end

