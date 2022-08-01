function [problem, state] = equationsWGVEbasicSens(model, state0, state, dt, drivingForces, varargin)
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
   opt = struct('reverseMode', false, ...
                'resOnly', false, ...
                'iteration', -1);

   opt = merge_options(opt, varargin{:});

   assert(isempty(drivingForces.src)); % unsupported

   s  = model.operators;
   f  = model.fluid;

   % Extract the current and previous values of all variables to solve for
   [p, sG, sGmax, wellSol,dz,rhofac,permfac,porofac] = model.getProps(state , 'pressure', 'sg', 'sGmax', 'wellsol', 'dz','rhofac','permfac','porofac');
   [p0, sG0,wellSol0]               = model.getProps(state0, 'pressure', 'sg','wellSol');

   [wellVars, wellVarNames, wellMap] = model.FacilityModel.getAllPrimaryVariables(wellSol);
   % Initialization of independent variables

   if ~opt.resOnly
      % ADI variables needed since we are not only computing residuals
      if ~opt.reverseMode
         [p, sG, wellVars{:}, dz, rhofac, permfac, porofac] = initVariablesADI(p, sG, wellVars{:}, dz,rhofac,permfac,porofac);
      else
         wellVars0 = model.FacilityModel.getAllPrimaryVariables(wellSol0);
         dzw = zeros(size(dz));
         dew = zeros(size(rhofac));
         [p0, sG0, wellVars0{:}, ~,~,~,~] = initVariablesADI(p0, sG0, wellVars0{:}, dzw,dew,dew,dew); %#ok
      end
   end
   primaryVars = {'pressure', 'sG', wellVarNames{:}};   
  
   sW  = 1 - sG;  % for ease of reading, we define an explicit variable
   sW0 = 1 - sG0; % also for water saturation

   % Preparing various necessary, intermediate values

   % multiplier for mobilities
   [pvMult, transMult, mobMult, pvMult0] = getMultipliers(f, p, p0);
   
   %
   %pvMult = porofac*pvMult;
   %pvMult0 = state0.porofac*pvMult0;
   %pvMult0 = porofac*pvMult0;
   %transMult = permfac*transMult;
   %transMult0 = state0.permfac*transMult0;
   %
   
   % relative permeability
   [krW, krG] = model.evaluateRelPerm({sW, sG}, p, 'sGmax', sGmax);
   krW = krW * mobMult;
   krG = krG * mobMult;

   % Transmissibility
   trans = s.T * transMult;

   % Gravity gradient per face, including the gravity component resulting
   % from caprock geometry
   %gdz = model.getGravityGradient();
   s  = model.operators;
   Gt = model.G;
   g  = model.gravity(3); % @@ requires theta=0
   gdz = g * s.Grad(Gt.cells.z +dz);
   
   % CO2 phase pressure
   pW  = p;
   pW0 = p0;
   pG  = p + f.pcWG(sG, p, 'sGmax', sGmax,'rhofac',rhofac);

   % Evaluate water and CO2 properties
   [vW, bW, mobW, rhoW, upcw] = ...
       getPhaseFluxAndProps_WGVEsens(model, pW, pG, krW, trans, gdz, 'W', 0, 0,rhofac);%,permfac,porofac);
   [vG, bG, mobG, rhoG, upcg] = ...
       getPhaseFluxAndProps_WGVEsens(model, pW, pG, krG, trans, gdz, 'G', 0, 0,rhofac);%,permfac,porofac);
  
   if model.outputFluxes
       state = model.storeFluxes(state, vW, [], vG);
   end
   
   bW0 = f.bW(pW0);
   bG0 = f.bG(pW0); % Yes, using water pressure also for gas here

   % Multiply upstream b-factors by interface fluxes to obtain fluxes at
   % standard conditions
   bWvW = s.faceUpstr(upcw, bW) .* vW;
   bGvG = s.faceUpstr(upcg, bG) .* vG;


   % Setting up brine and CO2 equations

   % Water (Brine)
   eqs{1} = porofac*((s.pv / dt) .* (pvMult .* bW .* sW - pvMult0 .* bW0 .* sW0)) + permfac*s.Div(bWvW);

   % Gas (CO2)
   eqs{2} = porofac*((s.pv / dt) .* (pvMult .* bG .* sG - pvMult0 .* bG0 .* sG0)) + permfac*s.Div(bGvG);
   types = {'cell' , 'cell'};
   names = {'water', 'gas'};  
   % Include influence of boundary conditions
   eqs = addFluxesFromSourcesAndBCSens(model, ...
           eqs, {pW, pG}, {rhoW, rhoG}, {mobW, mobG}, {bW, bG}, {sW, sG}, drivingForces,permfac, dz);

       
   dissolved={};
   [eqs, names, types, state.wellSol] = ...
       model.insertWellEquations(eqs, names, types, wellSol0, wellSol, ...
       wellVars, wellMap, p, ...
       {mobW, mobG}, {rhoW, rhoG}, dissolved, ...
       {}, dt, opt);
    eqs{6}=dz-drivingForces.dz;
    eqs{7}=rhofac-drivingForces.rhofac;
    eqs{8}=permfac-drivingForces.permfac;
    eqs{9}=porofac-drivingForces.porofac;
   % Setting up problem
   primaryVars = {primaryVars{:}, 'dz', 'rhofac','permfac','porofac'};
   types = {types{:}, 'scell','mult','mult','mult'};
   names = {names{:}, 'geometry','mrho','mperm','mporo'};


   problem = LinearizedProblem(eqs, types, names, primaryVars, state, dt);

end
