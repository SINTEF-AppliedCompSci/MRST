function [problem, state] = equationsWGVEdisgas(model, state0, state, dt, drivingForces, varargin)
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

   opt = struct('reverseMode', false, 'resOnly', false, 'iteration', -1);
   opt = merge_options(opt, varargin{:});

   assert(isempty(drivingForces.src)); % unsupported
   W  = drivingForces.W;
   s  = model.operators;
   Gt = model.G;
   f  = model.fluid;

   % Extract the current and previous values of all variables to solve for
   [p, sG, sGmax, rs, wellSol] = model.getProps(state , 'pressure', 'sg', 'sGmax', 'rs', 'wellsol');
   [p0, sG0, sGmax0, rs0, wellSol0]      = model.getProps(state0, 'pressure', 'sg', 'sGmax', 'rs', 'wellSol');
   
   [wellVars, wellVarNames, wellMap] = model.FacilityModel.getAllPrimaryVariables(wellSol);

   % Identifying cells with CO2-saturated brine (needed for correct switching of equations)
   rsSat  = f.rsSat(p);
   s_tol  = 0; %(f.dis_rate == 0) * sqrt(eps); % small nonzero if instantaneous dissolution model
   isSat  = (sG > s_tol) | (rs > rsSat);
   %isSat0 = (sG0 > s_tol);

   % Initialization of independent variables

   if ~opt.resOnly
      % ADI variables needed since we are not only computing residuals
      if ~opt.reverseMode
         [p, sG, sGmax, rs, wellVars{:}] = initVariablesADI(p, sG, sGmax, rs, wellVars{:});
      else
         wellVars0 = model.FacilityModel.getAllPrimaryVariables(wellSol0);
         [p0, sG0, sGmax0, rs0, wellVars0{:}] = initVariablesADI(p0, sG0, sGmax0, rs0, wellVars0{:});
      end
   end

   sW  = 1 - sG;  % for ease of reading, we define an explicit variable
   sW0 = 1 - sG0; % also for water saturation


   % Preparing various necessary, intermediate values

   % multiplier for mobilities
   [pvMult, transMult, mobMult, pvMult0] = getMultipliers(f, p, p0);

   % relative permeability
   [krW, krG] = model.evaluateRelPerm({sW, sG}, p, 'sGmax', sGmax);
   krW = krW * mobMult;
   krG = krG * mobMult;

   % Transmissibility
   trans = s.T .* transMult;

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

   if model.outputFluxes
       state = model.storeFluxes(state, vW, [], vG);
   end
   
   bW0 = f.bW(pW0);
   bG0 = f.bG(pW0); % Yes, using water pressure also for gas here

   % Multiply upstream b-factors by interface fluxes to obtain fluxes at
   % standard conditions
   bWvW   = s.faceUpstr(upcw, bW) .* vW;
   bGvG   = s.faceUpstr(upcg, bG) .* vG;
   rsbWvW = s.faceUpstr(upcw, rs) .* bWvW;

   % Setting up brine and CO2 equations
   eqs = cell(1, 4);
   % Water (Brine) conservation
   eqs{1} = (s.pv / dt) .* (pvMult .* bW .* sW - pvMult0 .* bW0 .* sW0) + s.Div(bWvW);

   % Gas (CO2) conservation
   eqs{2} = (s.pv / dt) .* (pvMult  .* (bG  .* sG  + rs  .* bW  .* sW) -   ...
                            pvMult0 .* (bG0 .* sG0 + rs0 .* bW0 .* sW0)) + ...
            s.Div(bGvG + rsbWvW);

   
   
    % Add in any fluxes / source terms prescribed as boundary conditions.
    dissolved = {{[], []}, {rs, []}};
    rho = {rhoW, rhoG};
    mob = {mobW, mobG};
    sat = {sW, sG};
    
   primaryVars = {'pressure' , 'sG'   , 'sGmax' , 'rs', wellVarNames{:}};
   types = {'cell'           , 'cell' , 'cell'  , 'cell'  };
   names = {'water'          , 'gas'  , 'dissol' , 'sGmax'};

    [eqs, state, src] = addBoundaryConditionsAndSources(model, eqs, names, types, state, ...
                                                     {pW, p}, sat, mob, rho, ...
                                                     dissolved, {}, ...
                                                     drivingForces);

    % Finally, add in and setup well equations
    [eqs, names, types, state.wellSol, qWell] = model.insertWellEquations(eqs, names, types, wellSol0, wellSol, wellVars, wellMap, p, mob, rho, dissolved, {}, dt, opt);
    wc = vertcat(drivingForces.W.cells);

   % Setting up dissolution equations
   eqs(3:4) = compute_dissolution_equations(model, Gt, f, sG, sG0, sGmax, ...
                                            sGmax0, bG, sW, sW0, bW, bW0, p, ...
                                            rhoW, mobW, rs, rs0, rsSat, isSat, ...
                                            pvMult, pvMult0, wc, qWell.phaseVolume{1}, s, ...
                                            rsbWvW, src, dt);

   % Rescaling non-well equations
   for i = 1:4
      eqs{i} = eqs{i} * dt / year;
   end


   problem = LinearizedProblem(eqs, types, names, primaryVars, state, dt);

end

% ============================================================================

function eqs = compute_dissolution_equations(model, Gt, f, sG, sG0, sGmax, ...
                                             sGmax0, bG, sW, sW0, bW, bW0, p, ...
                                             rhoW, mobW, rs, rs0, rsSat, isSat, ...
                                             pvMult, pvMult0, wc, cqsw, s, ...
                                             rsbWvW, src, dt)
   if f.dis_rate == 0

      % Instantaneous dissolution model
      eqs{1}        = sG + 0 * sG0;
      eqs{1}(isSat) = rs(isSat) - rsSat(isSat);
      eqs{2}        = sGmax - max(sGmax0, sG);

   else

      % Rate-driven dissolution model

      rate = f.dis_rate .* s.pv./Gt.cells.H; % rate per area multiplied
                                             % by CO2/brine interface area
                                             % in cell

      small = 3e-3;
      smooth_to_zero = @(x) value(x<small) .* ((x/small).^2 .* (- 2 * x / small + 3)) + value(x>=small) * 1;
      
      s_fac = smooth_to_zero(sG); % approximately one, but goes to 0 for very small values of sG
      rs_eps = (rsSat - rs) ./ rsSat;
      rs_fac = smooth_to_zero(rs_eps);
      rate = rate .* s_fac .* rs_fac; % smoothly turn down dissolution rate when
                                      % sG goes to 0, or when dissolved value
                                      % approaches maximum.

      % Conservation equation for dissolved CO2
      eqs{1} = (s.pv/dt) .* (pvMult  .* rs  .* bW  .* sW - pvMult0 .* rs0 .* bW0 .* sW0) + ...
               s.Div(rsbWvW) - rate;

      % accounting for dissolved gas in brine extracted from producer wells
      cqsw(cqsw > 0) = 0; % only take into account producer wells (inflow
                          % assumed not to contain any disgas)
      eqs{1}(wc) = eqs{1}(wc) - cqsw;

      % Taking boundary conditions and sources into account 
      
      flds = fieldnames(src);
      [vals, cells] = deal(cell(1, numel(flds)));
      for i = 1:numel(flds)
          obj = src.(flds{i});
          bcCells = obj.sourceCells;
          if ~isempty(bcCells)
              eqs{1}(bcCells) = eqs{1}(bcCells) - rs(bcCells).*obj.phaseMass{1}./model.fluid.rhoWS;
          end
          vals{i} = rs(bcCells).*obj.phaseMass{1}./model.fluid.rhoWS;
          cells{i} = bcCells;
      end

      % Modify conservation equation to ensure dissolved CO2 remains within
      % valid boundaries

      % Computing minimum allowed residual saturation per cell (the upper region where
      % both CO2 and brine phase is present should have saturated brine)
      min_rs = minRs(p, sG, sGmax, f, Gt);

      % Computing a hypothetic equation involving only minimum allowed
      % saturation.  We use this equation to identify where the real solution
      % would go below minimal allowed saturation.  Lock saturation of the
      % corresponding cells to the minimum value (but not lower).
      tmp = (s.pv / dt) .* ( value(pvMult) .* (value(min_rs) .* value(bW) ) - ...
                             pvMult0 .* (rs0 .* bW0 .* sW0 ) ) + ...
            s.Div(value(rsbWvW)) - value(rate);
      tmp(wc) = tmp(wc) - value(cqsw);
      for i = 1:numel(vals)
          tmp(cells{i}) = tmp(cells{i}) - value(vals{i});
      end

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
      is_sat = (value(rs) >= value(rsSat)) & (eqs{1} < 0);
      eqs{1}(is_sat) = (rs(is_sat) - rsSat(is_sat)) .* s.pv(is_sat)/dt;


      % Equation for changes to sGmax

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
             pvMult .* value(bG) .* (value(sG) - sGmax0) * ...
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
