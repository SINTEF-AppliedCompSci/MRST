function [problem, state] = equationsWGVEdisgas(model, state0, state, dt, drivingForces, varargin)

   opt = struct('reverseMode', false, 'resOnly', false, 'iteration', -1);
   opt = merge_options(opt, varargin{:});

   assert(isempty(drivingForces.src)); % unsupported
   W  = drivingForces.W;
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
   %isSat0 = (sG0 > s_tol);

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
   [vW, bW, mobW, rhoW, upcw] = ...
       getPhaseFluxAndProps_WGVE(model, pW, pG, krW, trans, gdz, 'W', rs, 0);
   [vG, bG, mobG, rhoG, upcg] = ...
       getPhaseFluxAndProps_WGVE(model, pW, pG, krG, trans, gdz, 'G', 0, 0);

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
   eqs = addFluxesFromSourcesAndBC(model,...
            eqs, {pW, pG}, {rhoW, rhoG}, {mobW, mobG}, {bW, bG}, {sW, sG}, drivingForces);
   eqs_rs{1} = eqs{1} * 0;
   eqs_rs{2} = eqs_rs{1};
   eqs_rs = addFluxesFromSourcesAndBC(model,...
            eqs_rs, {pW, 0*pG}, {rhoW, 0*rhoG}, {mobW, 0 * mobG}, {rs .* bW, ...
                       0*bG}, {sW, 0*sG}, drivingForces);
   eqs{2} = eqs{2} + eqs_rs{1};

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
         cqs = {0, 0}; % no in/outflow from any well
         wc  = [];     % no wellcells
         [eqs(5:7), names(5:7), types(5:7)] = ...
             wm.createReverseModeWellEquations(model, state0, wellSol, p0);%#ok
      end
   else
      eqs(5:7) = {bhp, bhp, bhp}; % empty ADIs
   end

   %% Setting up dissolution equations
   eqs(3:4) = compute_dissolution_equations(model, Gt, f, sG, sG0, sGmax, ...
                                            sGmax0, bG, sW, sW0, bW, bW0, p, ...
                                            rhoW, mobW, rs, rs0, rsSat, isSat, ...
                                            pvMult, pvMult0, wc, cqs{1}, s, ...
                                            rsbWvW, drivingForces, dt);

   %% Rescaling non-well equations
   for i = 1:4
      eqs{i} = eqs{i} * dt / year;
   end

   %% Setting up problem
   primaryVars = {'pressure' , 'sG'   , 'sGmax' , 'rs'    , 'qWs'      , 'qGs'          , 'bhp'       };
   types = {'cell'           , 'cell' , 'cell'  , 'cell'  , 'perf'     , 'well'         , 'perf'      };
   names = {'water'          , 'gas'  , 'dissol' , 'sGmax', 'gasWells' , 'closureWells' , 'waterWells'};
   if isempty(W)
      % Remove names/types associated with wells, as no well exist
      types = types(1:4);
      names = names(1:4);
   end

   problem = LinearizedProblem(eqs, types, names, primaryVars, state, dt);

end

% ============================================================================

function eqs = compute_dissolution_equations(model, Gt, f, sG, sG0, sGmax, ...
                                             sGmax0, bG, sW, sW0, bW, bW0, p, ...
                                             rhoW, mobW, rs, rs0, rsSat, isSat, ...
                                             pvMult, pvMult0, wc, cqsw, s, ...
                                             rsbWvW, drivingForces, dt)
   if f.dis_rate == 0

      %% Instantaneous dissolution model
      eqs{1}        = sG + 0 * sG0;
      eqs{1}(isSat) = rs(isSat) - rsSat(isSat);
      eqs{2}        = sGmax - max(sGmax0, sG);

   else

      %% Rate-driven dissolution model

      rate = f.dis_rate .* s.pv./Gt.cells.H; % rate per area multiplied
                                             % by CO2/brine interface area
                                             % in cell

      small = 3e-3;
      smooth_to_zero = @(x) double(x<small) .* ((x/small).^2 .* (- 2 * x / small + 3)) + double(x>=small) * 1;
      
      s_fac = smooth_to_zero(sG); % approximately one, but goes to 0 for very small values of sG
      rs_eps = (rsSat - rs) ./ rsSat;
      rs_fac = smooth_to_zero(rs_eps);
      rate = rate .* s_fac .* rs_fac; % smoothly turn down dissolution rate when
                                      % sG goes to 0, or when dissolved value
                                      % approaches maximum.

      %% Conservation equation for dissolved CO2
      eqs{1} = (s.pv/dt) .* (pvMult  .* rs  .* bW  .* sW - pvMult0 .* rs0 .* bW0 .* sW0) + ...
               s.Div(rsbWvW) - rate;

      % accounting for dissolved gas in brine extracted from producer wells
      cqsw(cqsw > 0) = 0; % only take into account producer wells (inflow
                          % assumed not to contain any disgas)
      eqs{1}(wc) = eqs{1}(wc) - cqsw;

      % Taking boundary conditions and sources into account 
      eqs_rs = {eqs{1} * 0, eqs{1} * 0};
      eqs_rs = addFluxesFromSourcesAndBC(model, eqs_rs, {p, 0*p}, {rhoW, 0*rhoW}, ...
                                         {mobW, 0*mobW}, {rs .* bW, 0*bW}, {sW, 0*sW}, ...
                                         drivingForces);
      eqs{1} = eqs{1} + eqs_rs{1};

      %% Modify conservation equation to ensure dissolved CO2 remains within
      % valid boundaries

      % Computing minimum allowed residual saturation per cell (the upper region where
      % both CO2 and brine phase is present should have saturated brine)
      min_rs = minRs(p, sG, sGmax, f, Gt);

      % Computing a hypothetic equation involving only minimum allowed
      % saturation.  We use this equation to identify where the real solution
      % would go below minimal allowed saturation.  Lock saturation of the
      % corresponding cells to the minimum value (but not lower).
      tmp = (s.pv / dt) .* ( double(pvMult) .* (double(min_rs) .* double(bW) ) - ...
                             pvMult0 .* (rs0 .* bW0 .* sW0 ) ) + ...
            s.Div(double(rsbWvW)) - double(rate);
      tmp(wc) = tmp(wc) - double(cqsw);
      tmp = tmp + double(eqs_rs{1});

      ilow = tmp > -sqrt(eps); % If so, then min_rs is larger than the
                               % solution of eqs{1}, in other words, the
                               % solution of eqs{1} is smaller than the
                               % allowed value.  We have to modify eqs{1}.
      if any(ilow)
         % force value of 'rs' in these cells to equal 'min_rs'
         eqs{1}(ilow) = rs(ilow) .* sW(ilow) - min_rs(ilow);
         eqs{1}(ilow) = eqs{1}(ilow) .* s.pv(ilow)/dt;
      end

      % Identify cells with fully-saturated brine, and ensure saturation does
      % not rise further
      is_sat = (double(rs) >= double(rsSat)) & (eqs{1} < 0);
      eqs{1}(is_sat) = (rs(is_sat) - rsSat(is_sat)) .* s.pv(is_sat)/dt;


      %% Equation for changes to sGmax

      % Default equation: force 'sGmax' to equal 'sG'
      eqs{2} = (sGmax - sG) .* (s.pv / dt);

      % In the following two special cases, sGmax will end up being higher
      % than sG, since any residual saturation is not fully eaten away by
      % dissolution (either due to rate being too slow, or due to saturation
      % of CO2 in brine being reached).
      tmp = (s.pv / dt) .* ...
            pvMult .* bG .* (sGmax - sGmax0) * ...
            f.res_gas ./ (1-f.res_water);
      tmp2 = (s.pv / dt) .* ...
             pvMult .* double(bG) .* (double(sG) - sGmax0) * ...
             f.res_gas ./ (1 - f.res_water) + rate;

      % Special case 1: New state remains UNSATURATED, but dissolution too slow to
      % deplete all residual saturation
      ix = (tmp2 < 0) & ~is_sat ;
      if any(ix)
         eqs{2}(ix) = tmp(ix) + rate(ix);
      end

      % Special case 2: New state reaches SATURATED value before all residual
      % saturation has been depleted
      ix = is_sat & (sGmax < sGmax0) & (sGmax > sG);
      if any(ix)
         eqs{2}(ix) = tmp(ix) + (rs(ix) - rs0(ix));
      end
   end
end
