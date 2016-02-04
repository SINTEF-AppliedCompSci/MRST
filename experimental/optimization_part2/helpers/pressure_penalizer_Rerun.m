function obj = pressure_penalizer_Rerun(model, states, schedule, penalty, plim, varargin)
% states.pressure is a scalar field (size of domain).
% schedule is only used for time steps.
% penalty is a scalar.
% plim is a scalar.


% format of objective function:
%   obj = max(0, sign(p - plim)) * penalty * (p - plim)^2
% obj is computed for each time step (for both inj and mig periods).


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
   maxp    = zeros(numSteps, 1);
   k       = 2;
   for step = 1:numSteps     
      state = states{tSteps(step)}; %@@ +1?      
      p = state.pressure;
      if opt.ComputePartials
        sG = state.s(:,2);   % place holders
        nW = numel(schedule.control(1).W);
        pBHP = zeros(nW, 1); % place holders
        qGs = pBHP;          % place holders
        qWs = pBHP;          % place holders
        [p, ~, ~, ~, ~] = initVariablesADI(p, sG, qWs, qGs, pBHP);
      end
      dt = dts(step); % how to scale obj value with dt?
      %maxp      = max(p); % a scalar for now
      tmp = - max(0, sign(p - 0.60*plim)) .* penalty .* (p - 0.60*plim).^k;
      tmp = tmp .* model.G.cells.volumes;
      obj{step} = sum( tmp )./sum(model.G.cells.volumes);
      %obj{step} = obj{step} * dt; % @@
      if (tSteps(step) == num_timesteps)
      % no need to compute another portion of the obj fun here
         if ~opt.ComputePartials
            %fprintf('Max pressure reached: %f (Pascals)\n', max(p) ); % @@ this is only the max pressure reached at last time step!
            %fprintf('Proximity to PLimit: %f (percent)\n', maxp/(0.75*plim) * 100 );
            %fprintf('Score: %f ([Cp*Pascals])\n\n', max(0, sign(maxp - 0.75*plim)) * penalty * (maxp - 0.75*plim).^k );
          end
      end
      %obj{step} = obj{step};  %@@ how to scale this value to match with other obj value?
   end
end
