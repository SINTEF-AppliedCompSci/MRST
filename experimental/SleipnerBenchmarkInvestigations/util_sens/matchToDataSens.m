function obj = matchToDataSens(model, wellSols, states, schedule, plumes, varargin)
% compute match between 'simulated' and 'observed' CO2 heights

% obj is the L2-norm of the difference between 'observed' CO2 heights
% (plumes.h) and 'simulated' CO2 heights

   
   opt.ComputePartials = false;
   opt.tStep = [];
   opt.pl = [];
   opt = merge_options(opt, varargin{:});
  
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
   for step = 1:numSteps
      sol = wellSols{tSteps(step)};
      state = states{tSteps(step)}; %@@ +1?
      nW = numel(sol);
      pBHP = zeros(nW, 1); % place holder
      qGs = vertcat(sol.qGs);
      qWs = vertcat(sol.qWs);
      p = state.pressure;
      sG = state.s(:,2);
      dz = state.dz;
      rhofac=state.rhofac;
      permfac=state.permfac;
      porofac=state.porofac;
      if opt.ComputePartials
         [p, sG, pBHP, qWs, qGs, dz, rhofac, permfac, porofac] = initVariablesADI(p, sG, pBHP, qWs, qGs, dz, rhofac, permfac, porofac);%#ok
      end
      %dt = dts(step);
      %injInx = (vertcat(sol.sign) > 0);
      %obj{step} = dt * spones(ones(1, nW)) * ((1-penalty) * injInx .* qGs);
      
      if(~isempty(plumes{tSteps(step)}))
          obj{step} = matchToCo2Surface(sG,plumes{tSteps(step)},model.G,model.fluid);
          obj{step} = obj{step} + sum(dz.^2.*model.G.cells.volumes); 
      else
         obj{step} = 0*qGs;
      end      
   end
end

