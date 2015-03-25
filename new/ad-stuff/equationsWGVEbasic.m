function [problem, state] = equationsWGVEbasic(state0, state, model, dt, drivingForces, varargin)
   
   opt = struct('reverseMode', false, ...
                'resOnly', false);
   
   assert(isempty(drivingForces.src)); % unsupported
   W  = drivingForces.Wells;
   bc = drivingForces.bc; 
   s  = model.operators; 
   Gt = model.G; 
   f  = model.fluid; 
   
   % Extract the current and previous values of all variables to solve for
   [p, sG, sGmax, wellSol] = model.getProps(state , 'pressure', 'sg', 'sGmax', 'wellsol');
   [p0, sG0]               = model.getProps(state0, 'pressure', 'sg');
   
   % Stack well-related variables of the same type together
   bhp = vertcat(wellSol.bhp);
   qWs = vertcat(wellSol.qWs);
   qGs = vertcat(wellSol.qGs);
   
   %% Initialization of independent variables
   
   if ~opt.resOnly
      % ADI variables needed since we are not only computing residuals
      if ~opt.reverseMode
         [p, sG, bhp, qWs, qGs] = initVariablesADI(p, sG, bhp, qWs, qGs);
      else
         zw = zeros(size(bhp)); % dummy
         [p0, sG0, ~, ~, ~] = initVariablesADI(p0, sG0, zw, zw, zw);
      end
   end
   
   sW  = 1 - sG;  % for ease of reading, we define an explicit variable 
   sW0 = 1 - sG0; % also for water saturation 
   
   %% Preparing various necessary, intermediate values
   
   % multiplier for mobilities
   [pvMult, transMult, mobMult, pvMult0] = getMultipliers(f, p, p0);
   
   % relative permeability
   [krW, krG] = model.evaluteRelPerm({sW, sG}, p, 'sGmax', sGmax);
   krW = krW * mobMult;
   krG = krG * mobMult;
   
   % Transmissibility
   trans = s.T * transMult;
   
   % Gravity gradient per face, including the gravity component resulting
   % from caprock geometry  
   gdz = model.getGravityGradient();

   % CO2 phase pressure
   pg = p + pcWG(sG, p, 'sGmax', sGmax);
   
   % Evaluate water and CO2 properties 
   [vW, bW, mobW, rhoW, upcw, dpW] = getPhaseFluxAndProps_WGVE(model, pw, pg, krW, trans, gdz, 'W');
   [vG, bG, mobG, rhoG, upcg, dpG] = getPhaseFluxAndProps_WGVE(model, pw, pg, krG, trans, gdz, 'G');
   bW0 = f.bW(pW);
   bG0 = f.bG(pW); % Yes, using water pressure also for gas here
   
   % Multiply upstream b-factors by interface fluxes to obtain fluxes at
   % standard conditions
   bWvW = s.faceUpstr(upcw, bW) .* vW;
   bGvG = s.faceUpstr(upcg, bG) .* vG;
   
   
   %% Setting up equations 
   
   % Gas (CO2)
   eqs{1} = (s.pv / dt) .* (pvMult .* bW .* sW - pvMult0 .* bW0 .* sW0) + s.div(bWvW);
   
   % Water (Brine)
   eqs{2} = (s.pv / dt) .* (pvMult .* bG .* sG - pvMult0 .* bG0 .* sG0) + s.div(bGvG);
   
   
   % ADD BOUNDARY STUFF
   
   % SETUP WELL EQUATIONS
   
   
   %% Setting up problem
   types = {'cell', 'cell', 'cell', well};
   names = {'water', 'gas', wellstuff};
   primaryVars = {'pressure', 'sG', 'bhp', 'qWs', 'qGs'};
   problem = LinearizedProblem(eqs, types, names, primaryVars, state, dt);
   
end
