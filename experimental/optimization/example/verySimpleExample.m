function [Gt, optim] = verySimpleExample()

   moduleCheck('ad-core');
   gravity on;
   
   ref_pressure = 100 * barsa;
   rate = 1; %0.18;%5;%sqrt(eps);
   
   res = 33;
   G = computeGeometry(cartGrid([res res 1], [1 1 0.02] * kilo * meter));
   Gt = topSurfaceGrid(G);
   rock.perm = 100 * milli * darcy * ones(Gt.cells.num, 1);
   rock.poro = 0.2 * ones(Gt.cells.num, 1);
   rock = averageRock(rock, Gt);
   
   fluid = makeVEFluid(Gt, rock, 'simple', ...
                       'fixedT', 40 + 275.25, ...
                       'co2_rho_pvt', [], ...
                       'wat_rho_pvt', [], ...
                       'pvMult_fac',  0);
                       
   wcell = (res^2 + 1)/2; %5;
   W = addWell([], Gt, rock, wcell    , ...
               'Type'   , 'rate'      , ...
               'Val'    , rate        , ...
               'Radius' , 0.3         , ...
               'Comp_i' , [0, 1]      , ...
               'name'   , ['I']);
   W = addWell(W, Gt, rock, 1    , ...
               'Type'   , 'rate'      , ...
               'Val'    , rate        , ...
               'Radius' , 0.3         , ...
               'Comp_i' , [0, 1]      , ...
               'name'   , ['I']);
   
   schedule.control.W = W;
   bfaces = identifyBoundaryFaces(Gt);
   schedule.control.bc = addBC([], bfaces, 'pressure', ref_pressure, ...
                                   'sat', [1 0]);
   schedule.step.val = ones(1, 10) * day;
   schedule.step.control = ones(1, 10);
                                   
   initState.pressure = ref_pressure * ones(Gt.cells.num, 1);
   initState.s = repmat([1 0], Gt.cells.num, 1);
   initState.sGmax = initState.s(:,2);
   
   model = CO2VEBlackOilTypeModel(Gt, rock, fluid);
   min_rates = 0.01 * ones(numel(W), 1);
   max_rates = 30   * ones(numel(W), 1);
   
   %obj_wrapper = @(ws, st, sch, varargin) one_step(model, ws, st, sch, varargin{:});
   %obj_wrapper = @(ws, st, sch, varargin) simple_leak(model, ws, st, sch, varargin{:});
   %obj_wrapper = @(ws, st, sch, varargin) try_again(model, ws, st, sch, varargin{:});   
   obj_wrapper = @(ws, st, sch, varargin) nico(model, ws, st, sch, varargin{:});   
   
   [optim, init, history] = optimizeRates(initState, model, schedule, min_rates, ...
                                          max_rates, 'obj_fun', obj_wrapper);
   % [optim, init, history] = optimizeRates(initState, model, schedule, min_rates, ...
   %                                        max_rates);
   
end

function obj = nico(model, wellSols, states, schedule, varargin)
   
   opt.ComputePartials = false;
   opt.tStep = [];
   opt.penalty = 10;
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
   krull = 0;
   for step = 1:numSteps
      sol = wellSols{tSteps(step)};
      state = states{tSteps(step)}; %@@ +1?
      nW = numel(sol);
      pBHP = zeros(nW, 1); % place holder
      qGs = vertcat(sol.qGs);
      qWs = vertcat(sol.qWs);
      p = state.pressure;
      sG = state.s(:,2);
      if opt.ComputePartials
         [p, sG, pBHP, qWs, qGs] = initVariablesADI(p, sG, pBHP, qWs, qGs);
      end
      dt = dts(step);
      injInx = (vertcat(sol.sign) > 0);
      obj{step} = dt * spones(ones(1, nW)) * ((1-opt.penalty) * injInx .* qGs);

      krull = krull +  dt * spones(ones(1, nW)) * ( injInx .* qGs);
      
      if (tSteps(step) == num_timesteps)
         bG = model.fluid.bG(p);
         pvol = model.G.cells.volumes .* model.G.cells.H .* model.rock.poro;      
         vol = ones(1, model.G.cells.num) * (pvol .* model.fluid.pvMultR(p) .* bG .* sG);
         obj{step} = obj{step} + opt.penalty * vol;
         if ~opt.ComputePartials
            fprintf('Total injected: %f\n', double(krull));
            fprintf('Total leaked: %f\n', double(krull - vol));
            fprintf('Score: %f\n\n', double(krull) - opt.penalty * double(krull-vol));
         end
      end
      obj{step} = obj{step} * model.fluid.rhoGS/1e12;
   end
end


function obj = try_again(model, wellSols, states, schedule, varargin)

   opt.ComputePartials = false;
   opt.tStep = [];
   opt.penalty = 5;
   opt = merge_options(opt, varargin{:});
   
   pvol = model.G.cells.volumes .* model.G.cells.H .* model.rock.poro;
   if isfield(model.fluid, 'pvMultR')
      pvMult = model.fluid.pvMultR;
   else
      pvMult = @(p) p*0+1;
   end
   numSteps = numel(states);
   injected = zeros(numSteps, 1);
   in_place = zeros(numSteps, 1);
   leaked   = zeros(numSteps, 1);
   
   for i = 1:numSteps
      
      st = states{i};
      dt = schedule.step.val(i);
      P  = st.pressure;
      sG = st.s(:,2);
      qGs = [wellSols{i}.val];
      
      injected(i) = sum(qGs * dt);
      in_place(i) = sum(pvol .* pvMult(P) .* sG .* model.fluid.bG(P));
      if i == 1
         prev_in_place = 0;
      else
         prev_in_place = in_place(i-1);
      end
      
      leaked(i) = injected(i) - (in_place(i) - prev_in_place);
   end
   
   if isempty(opt.tStep)
      assert(~opt.ComputePartials)
      obj = cell(numSteps, 1);
      for i = 1:numSteps
         obj{i} = injected(i) - opt.penalty * leaked(i);
      end
   else
      assert(opt.ComputePartials)
      tStep = opt.tStep;
      dt = schedule.step.val(tStep);
      st = states{tStep};
      P  = st.pressure;
      sG = st.s(:,2);
      qGs = [wellSols{tStep}.val];

      dummy = zeros(size(qGs));
      [P, sG, ~, ~, qGs] = initVariablesADI(P, sG, dummy, dummy, qGs);
      
      if numel(double(qGs)) > 1
         cur_injected = sum(qGs * dt);
      else
         cur_injected = qGs * dt;
      end
      cur_in_place =  sum(pvol .* pvMult(P) .* sG .* model.fluid.bG(P));
      if tStep == 1
         cur_leaked = cur_injected - cur_in_place;
      else 
         cur_leaked = cur_injected - (cur_in_place - in_place(tStep - 1));
      end
      
      obj = cur_injected - opt.penalty * cur_leaked;
   
   end   
obj
end


function obj = one_step(model, wellSols, states, schedule, varargin)
   
   opt.ComputePartials = false;
   opt.tStep = [];
   opt.penalty = 10;
   opt = merge_options(opt, varargin{:});
   
   pvol = model.G.cells.volumes .* model.G.cells.H .* model.rock.poro;
   if isfield(model.fluid, 'pvMultR')
      pvMult = model.fluid.pvMultR;
   else
      pvMult = @(p) p*0+1;
   end
   numSteps = numel(states);
   assert(numSteps == 1);
   dt = schedule.step.val;
   
   P = states{1}.pressure;
   sG = states{1}.s(:,2);

   qGs = [wellSols{1}.val];
   
   dummy = zeros(size(qGs));
   [P, sG, ~, ~, qGs] = initVariablesADI(P, sG, dummy, dummy, qGs);
   
   assert(numel(double(qGs))== 1);
   
   injected = qGs * dt;
   
   remaining = sum(pvol .* pvMult(P) .* sG .* model.fluid.bG(P));       
   
   leaked = injected - remaining;
   if ~opt.ComputePartials
      injected = double(injected);
      leaked = double(leaked);
   end
   
   obj = {injected - 10 * leaked};
end

   

function obj = simple_leak(model, wellSols, states, schedule, varargin)

   opt.ComputePartials = false;
   opt.tStep = [];
   opt.penalty = 10;
   opt = merge_options(opt, varargin{:});
   
   pvol = model.G.cells.volumes .* model.G.cells.H .* model.rock.poro;
   if isfield(model.fluid, 'pvMultR')
      pvMult = model.fluid.pvMultR;
   else
      pvMult = @(p) ones(size(p));
   end

   numSteps = numel(states);
   obj = cell(numSteps, 1);
   
   tot_injected = 0;
   for i = 1:numSteps
      obj{i} = 0;
      dt = schedule.step.val(i);
      tot_injected = tot_injected + sum([wellSols{i}.val]) * dt;
   end

   % Computing amount remaining at end
   fstate = states{end};
   Pfinal = fstate.pressure;
   sGfinal = fstate.s(:,2);

   tot_remaining = sum(pvol .* pvMult(Pfinal) .* sGfinal .* model.fluid.bG(Pfinal)); 
   
   leaked = tot_injected - tot_remaining;

   
   if ~isempty(opt.tStep)
      assert(opt.ComputePartials);
      st = states{opt.tStep};
      ws = wellSols{opt.tStep};
      
      P = st.pressure;
      sG = st.s(:,2);
      qGs = [ws.val];
      dt = schedule.step.val(opt.tStep);

      dummy = zeros(size(qGs));
      [P, sG, ~, ~, qGs] = initVariablesADI(P, sG, dummy, dummy, qGs);
      
      if numel(qGs.val) == 1
         cur_inj = qGs * dt;
      else
         cur_inj = sum(qGs) * dt;
      end
      
      % Adding relevant partial derivatives to the total injected amount and 
      tot_injected = tot_injected - double(cur_inj) + cur_inj;
      
      if opt.tStep == numSteps
         % giving 'tot_remaining' appropriate ADI derivatives
         tot_remaining = sum(pvol .* pvMult(P) .* sG .* model.fluid.bG(P));
      end
      
      % the 'leaked' variable will be recomputed to give it the right set of
      % derivatives 
      leaked  = tot_injected - tot_remaining;
      assert(leaked > -sqrt(eps));
   end
   
   obj{end} = tot_injected - 10 * leaked;%@@tot_injected - opt.penalty * leaked;   
   
   if ~isempty(opt.tStep)
      obj = obj{opt.tStep};
      if (opt.tStep ~= numSteps)
         obj = 0 * P(1); % let it have an ADI structure, although it is zero.
      end
   end
   
end




function obj = simple_obj(model, wellSols, states, schedule, varargin)
   
   opt.ComputePartials = false;
   opt.tStep = [];
   opt.penalty = 10;
   opt = merge_options(opt, varargin{:});
   
   pvol = model.G.cells.volumes .* model.G.cells.H .* model.rock.poro;
   if isfield(model.fluid, 'pvMultR')
      pvMult = model.fluid.pvMultR;
   else
      pvMult = @(p) ones(size(p));
   end
   
   
   numSteps = numel(states);
   obj = cell(numSteps, 1);
   for i = 1:numSteps
      dt = schedule.step.val(i);
      obj{i} = sum([wellSols{i}.val]) * dt;
   end

   total_injected = sum([obj{:}]);
   fstate = states{end};
   Pfinal = fstate.pressure;
   sGfinal = fstate.s(:,2);
   total_remaining = sum(pvol .* pvMult(Pfinal) .* sGfinal .* model.fluid.bG(Pfinal));
   leaked = total_injected - total_remaining;
   
   fprintf('injected: %f     leaked: %f\n', total_injected, leaked);
   % Adding penalty term to last timestep
   obj{end} = obj{end} - leaked * opt.penalty;
   
   fprintf('obj_end: %f\n\n', obj{end});

   if ~isempty(opt.tStep)
      assert(opt.ComputePartials)
      
      qGs = [wellSols{opt.tStep}.val];
      P   = states{opt.tStep}.pressure;
      sG  = states{opt.tStep}.s(:,2);
      % dummy vars
      %dummy_big = zeros(size(states{1}.pressure));
      dummy_small = zeros(size(qGs));
      [P, sG, ~, ~, qGs] = initVariablesADI(P, sG, dummy_small, dummy_small, qGs);
      obj = qGs;

      if (opt.tStep == numel(states))
         % This is the last timestep - we should also include the penalty
         % term
         total_injected = total_injected - double(qGs) + qGs; % making it ADI
         total_remaining = sum(pvol .* pvMult(P) .* sG .* model.fluid.bG(P));
         leaked = total_injected - total_remaining;
         
         obj = obj - leaked * opt.penalty;
      end
   end
   
   
end

   
function bfaces = identifyBoundaryFaces(Gt)
   
    bfaces = find(any(Gt.faces.neighbors==0,2)); 

end
