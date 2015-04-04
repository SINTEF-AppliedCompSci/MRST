function [problem, state] = equationsWGVEdisgas(model, state0, state, dt, drivingForces, varargin)

   opt = struct('reverseMode', false, 'resOnly', false, 'iteration', -1);
   opt = merge_options(opt, varargin{:});
   
   assert(isempty(drivingForces.src)); % unsupported
   W  = drivingForces.Wells;
   bc = drivingForces.bc; 
   s  = model.operators; 
   Gt = model.G; 
   f  = model.fluid; 
   
   % Extract the current and previous values of all variables to solve for
   [p, sG, sGmax, rs, wellSol] = model.getProps(state , 'pressure', 'sg', 'sGmax', 'rs', 'wellsol');
   [p0, sG0, sGmax0, rs0]      = model.getProps(state0, 'pressure', 'sg', 'sGmax', 'rs');
   
   % Stack well-related variables of the same type together
   bhp = vertcat(wellSol.bhp);
   qWs = vertcat(wellSol.qWs);
   qGs = vertcat(wellSol.qGs);   
   
   % Identifying cells with CO2-saturated brine (needed for correct switching of equations)
   rsSat  = f.rsSat(p);
   s_tol  = 0; %(f.dis_rate == 0) * sqrt(eps); % small nonzero if instantaneous dissolution model
   isSat  = (sG > s_tol) | (rs > rsSat);
   isSat0 = (sG0 > s_tol);
         
   %% Initialization of independent variables
   
   if ~opt.resOnly
      % ADI variables needed since we are not only computing residuals
      if ~opt.reverseMode
         [p, sG, sGmax, rs, qWs, qGs, bhp] = initVariablesADI(p, sG, sGmax, rs, qWs, qGs, bhp);
      else
         zw = zeros(size(bhp)); % dummy
         [p0, sG0, sGmax0, rs0, ~, ~, ~] = initVariablesADI(p0, sG0, sGmax0, rs0, zw, zw, zw);
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
   pW  = p;
   pW0 = p0;
   pG  = p + f.pcWG(sG, p, 'sGmax', sGmax);
   
   % Evaluate water and CO2 properties 
   [vW, bW, mobW, rhoW, upcw, dpW] = ...
       getPhaseFluxAndProps_WGVE(model, pW, pG, krW, trans, gdz, 'W', rs);
   [vG, bG, mobG, rhoG, upcg, dpG] = ...
       getPhaseFluxAndProps_WGVE(model, pW, pG, krG, trans, gdz, 'G', 0);

   bW0 = f.bW(pW0);
   bG0 = f.bG(pW0); % Yes, using water pressure also for gas here
   
   % Multiply upstream b-factors by interface fluxes to obtain fluxes at
   % standard conditions
   bWvW   = s.faceUpstr(upcw, bW) .* vW;
   bGvG   = s.faceUpstr(upcg, bG) .* vG;
   rsbWvW = s.faceUpstr(upcw, rs) .* bWvW;
   
   %% Setting up brine and CO2 equations 
   
   % Water (Brine) conservation
   eqs{1} = (s.pv / dt) .* (pvMult .* bW .* sW - pvMult0 .* bW0 .* sW0) + s.Div(bWvW);
   
   % Gas (CO2) conservation
   eqs{2} = (s.pv / dt) .* (pvMult  .* (bG  .* sG  + rs  .* bW  .* sW) -   ...
                            pvMult0 .* (bG0 .* sG0 + rs0 .* bW0 .* sW0)) + ...
            s.Div(bGvG + rsbWvW);
   
   % Include influence of boundary conditions (first line adds the phase
   % fluxes, second line adds the transport of gas in the water flux)
   eqs = addFluxesFromSourcesAndBC(model, ...
           eqs, {pW, pG}, {rhoW, rhoG}, {mobW, mobG}, {bW, bG}, {sW, sG}, drivingForces);
   eqs = addFluxesFromSourcesAndBC(model, ...
           eqs, {0, pW}, {0, rhoW}, {0, mobW}, {0, rs .* bW}, {0, sW}, drivingForces);
           
   %% Setting up dissolution equations
   eqs(3:4) = compute_dissolution_equations(f, sG, sG0, sGmax, sGmax0, rs, rsSat, isSat);
   
   %% Setting up well equations
   if ~isempty(W)
      wm = model.wellmodel;
      if ~opt.reverseMode
         wc = vertcat(W.cells);
         [cqs, weqs, ctrleqs, wc, state.wellSol] =                       ...
             wm.computeWellFlux(model, W, wellSol, bhp                 , ...
                                {qWs, qGs}                             , ...
                                p(wc)                                  , ...
                                [model.fluid.rhoWS, model.fluid.rhoGS] , ...
                                {bW(wc), bG(wc)}                       , ...
                                {mobW(wc), mobG(wc)}                   , ...
                                {sW(wc), sG(wc)}                       , ...
                                {}                                     , ...
                                'allowControlSwitching', false         , ... 
                                'nonlinearIteration', opt.iteration);
         % Store the separate well equations (relate well bottom hole
         % pressures to influx)
         eqs(5:6) = weqs;
         
         % Store the control equations (ensuring that each well has values
         % corresponding to the prescribed value)
          eqs{7} = ctrleqs;
         
         % Add source term to equations.  Negative sign may be surprising if
         % one is used to source terms on the right hand side, but this is
         % the equations on residual form
         eqs{1}(wc) = eqs{1}(wc) - cqs{1};
         eqs{2}(wc) = eqs{2}(wc) - cqs{2};
      
      else 
         [eqs(5:7), names(5:7), types(5:7)] = ...
             wm.createReverseModeWellEquations(model, state0, wellSol, p0);
      end
   end
   
   %% Rescaling non-well equations
   for i = 1:4
      eqs{i} = eqs{i} * dt / year;
   end
   
   
   %% Setting up problem
   primaryVars = {'pressure' , 'sG'   , 'sGmax' , 'rs'    , 'qWs'      , 'qGs'          , 'bhp'       };
   types = {'cell'           , 'cell' , 'cell'  , 'cell'  , 'perf'     , 'well'         , 'perf'      };
   names = {'water'          , 'gas'  , 'sGmax' , 'dissol', 'gasWells' , 'closureWells' , 'waterWells'};
   if isempty(W)
      % Remove names/types associated with wells, as no well exist
      types = types(1:4);
      names = names(1:4);
   end
      
   problem = LinearizedProblem(eqs, types, names, primaryVars, state, dt);
   
end

% ============================================================================

function eqs = compute_dissolution_equations(fluid, sG, sG0, sGmax, sGmax0, rs, rsSat, isSat)
   
   if fluid.dis_rate == 0
      % Instantaneous dissolution model
      eqs{1}        = sG + 0 * sG0;
      eqs{1}(isSat) = rs(isSat) - rsSat(isSat);
      eqs{2}        = sGmax - max(sGmax0, sG);
   else
      % Rate-driven dissolution model
      error('Rate-driven dissolution model currently unimplemented.');
   end
   
end
