function [problem, state] = equationsWGVEbasic(model, state0, state, dt, drivingForces, varargin)
%{
Copyright 2009-2016 SINTEF ICT, Applied Mathematics.

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

   opt = struct('reverseMode', false, ...
                'resOnly', false, ...
                'iteration', -1,...
                'adjointForm',false);

   opt = merge_options(opt, varargin{:});

   assert(isempty(drivingForces.src)); % unsupported
   W  = drivingForces.W;
   s  = model.operators;
   f  = model.fluid;

   % Extract the current and previous values of all variables to solve for
   if(opt.adjointForm)
      % use sGmax as primary variable
      [p, sG, sGmax, wellSol] = model.getProps(state , 'pressure', 'sg', 'sGmax', 'wellsol');
   else
      [p, sG, wellSol] = model.getProps(state , 'pressure', 'sg', 'wellsol');
   end
   [p0, sG0, sGmax0]               = model.getProps(state0, 'pressure', 'sg','sGmax');

   % Stack well-related variables of the same type together
   bhp = vertcat(wellSol.bhp);
   qWs = vertcat(wellSol.qWs);
   qGs = vertcat(wellSol.qGs);

   %% Initialization of independent variables

   if ~opt.resOnly
      % ADI variables needed since we are not only computing residuals
      if ~opt.reverseMode
         if(opt.adjointForm)
             [p, sG, sGmax, qWs, qGs, bhp] = initVariablesADI(p, sG, sGmax, qWs, qGs, bhp);
         else
            [p, sG, qWs, qGs, bhp] = initVariablesADI(p, sG, qWs, qGs, bhp);
            sGmax = max(sG, sGmax0);
         end
      else
          zw = zeros(size(bhp));
         if(opt.adjointForm)
            [p0, sG0, sGmax0, ~, ~, ~] = initVariablesADI(p0, sG0, sGmax0, zw, zw, zw);
         else
             warning('using non adjoint form for backward simulation: Hysteresis may be not properly handled')
            [p0, sG0, ~, ~, ~] = initVariablesADI(p0, sG0, zw, zw, zw);
            sGmax = max(sG, sGmax0);
         end
      end
   elseif ~opt.adjointForm
      % sGmax was not yet defined.  We are probably in the process of
      % evaluating residuals after non-convergence.  Set it here to avoid
      % problems further down
      sGmax = max(sG, sGmax0);
   end
   

   sW  = 1 - sG;  % for ease of reading, we define an explicit variable
   sW0 = 1 - sG0; % also for water saturation

   %% Preparing various necessary, intermediate values

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
       getPhaseFluxAndProps_WGVE(model, pW, pG, krW, trans, gdz, 'W', 0, 0);
   [vG, bG, mobG, rhoG, upcg] = ...
       getPhaseFluxAndProps_WGVE(model, pW, pG, krG, trans, gdz, 'G', 0, 0);
   bW0 = f.bW(pW0);
   bG0 = f.bG(pW0); % Yes, using water pressure also for gas here

   % Multiply upstream b-factors by interface fluxes to obtain fluxes at
   % standard conditions
   bWvW = s.faceUpstr(upcw, bW) .* vW;
   bGvG = s.faceUpstr(upcg, bG) .* vG;


   %% Setting up brine and CO2 equations

   % Water (Brine)
   eqs{1} = (s.pv / dt) .* (pvMult .* bW .* sW - pvMult0 .* bW0 .* sW0) + s.Div(bWvW);

   % Gas (CO2)
   eqs{2} = (s.pv / dt) .* (pvMult .* bG .* sG - pvMult0 .* bG0 .* sG0) + s.Div(bGvG);

   % Include influence of boundary conditions
   eqs = addFluxesFromSourcesAndBC(model, ...
           eqs, {pW, pG}, {rhoW, rhoG}, {mobW, mobG}, {bW, bG}, {sW, sG}, drivingForces);

   if(opt.adjointForm)
       eqs{3}=sGmax-max(sG,sGmax0);
       eqshift=1;
   else
       eqshift=0;
   end
       
       

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
         eqs([3:4]+eqshift) = weqs;%#ok

         % Store the control equations (ensuring that each well has values
         % corresponding to the prescribed value)
         eqs{5+eqshift} = ctrleqs;

         % Add source term to equations.  Negative sign may be surprising if
         % one is used to source terms on the right hand side, but this is
         % the equations on residual form
         eqs{1}(wc) = eqs{1}(wc) - cqs{1};
         eqs{2}(wc) = eqs{2}(wc) - cqs{2};

      else
         [eqs([3:5]+eqshift), names([3:5]+eqshift), types([3:5]+eqshift)] = ...
             wm.createReverseModeWellEquations(model, state0.wellSol, p0);%#ok
      end
   else
      eqs([3:5]+eqshift) = {bhp, bhp, bhp}; %#ok empty ADIs
   end


   %% Setting up problem
   % original
   %primaryVars = {'pressure' , 'sG'  , 'qWs'      , 'qGs', 'bhp'};
   %types = {'cell'           , 'cell' , 'perf'       , 'perf'     , 'well'};
   %names = {'water'          , 'gas'  , 'waterWells' , 'gasWells' , 'closureWells'};
   %
   primaryVars = {'pressure' , 'sG'}; 
   types = {'cell'           , 'cell'};
   names = {'water'          , 'gas'};  
   %add hysteres variable,equation and name
   if(opt.adjointForm)
      primaryVars = {primaryVars{:}, 'sGmax'}; %#ok
      types = {types{:} , 'cell' };  %#ok   
      names = {names{:} , 'gasMax'}; %#ok
   end
   % add names of well equations
   primaryVars = {primaryVars{:}  , 'qWs'      , 'qGs', 'bhp'}; %#ok
   types = {types{:} , 'perf'       , 'perf'     , 'well'}; %#ok
   names = {names{:} , 'waterWells' , 'gasWells' , 'closureWells'}; %#ok
   
   if isempty(W)
      % Remove names/types associated with wells, as no well exist
      types = types(1:(2+eqshift));
      names = names(1:(2+eqshift));
   end

   problem = LinearizedProblem(eqs, types, names, primaryVars, state, dt);

end
