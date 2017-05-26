function [problem, state] = equationsOilWaterPolymer(state0, state, model, dt, ...
                                                     drivingForces, varargin)
%
% SYNOPSIS:
%   function [problem, state] = equationsOilWaterPolymer(state0, state, model, dt, drivingForces, varargin)
%
% DESCRIPTION: Assemble the linearized equations for an oil-water-polymer
% system, computing both the residuals and the Jacobians. Returns the result as
% an instance of the class LinearizedProblem which can be solved using instances
% of LinearSolverAD.
%
% A description of the modeling equations can be found in the directory
% ad-eor/docs.
%
%
% PARAMETERS:
%   state0        - State at previous times-step
%   state         - State at current time-step
%   model         - Model instance
%   dt            - time-step
%   drivingForces - Driving forces (boundary conditions, wells, ...)
%   varargin      - optional parameters
%
% RETURNS:
%   problem - Instance of LinearizedProblem
%   state   - Updated state variable (fluxes, mobilities and more can be
%             stored, the wellSol structure is also updated in case of control switching)
%
% EXAMPLE:
%
% SEE ALSO: LinearizedProblem, LinearSolverAD, equationsOilWater, OilWaterPolymerModel
%

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


% Get linearized problem for oil/water/polymer system with black oil-style
% properties
opt = struct('Verbose', mrstVerbose, ...
             'reverseMode', false,...
             'resOnly', false,...
             'iteration', -1);

opt = merge_options(opt, varargin{:});

W = drivingForces.W;
op = model.operators;

% Properties at current timestep
[p, sW, c, cmax, wellSol] = model.getProps(state, 'pressure', 'water', ...
    'polymer', 'polymermax', 'wellSol');

% Properties at previous timestep
[p0, sW0, c0, cmax0, wellSol0] = model.getProps(state0, 'pressure', 'water', ...
   'polymer', 'polymermax', 'wellSol');

[qWell, pBH, wellVars, wellVarNames, wellMap] = model.FacilityModel.getAllPrimaryVariables(wellSol);
% Initialize independent variables.
if ~opt.resOnly,
    % ADI variables needed since we are not only computing residuals.
    if ~opt.reverseMode,
        [p, sW, c, qWell{:}, pBH, wellVars{:}] = ...
            initVariablesADI(p, sW, c, qWell{:}, pBH, wellVars{:});
    else
        zw = zeros(size(pBH));
        [p0, sW0, c0, zw, zw, zw, zw] = ...
            initVariablesADI(p0, sW0, c0, zw, zw, zw, zw); %#ok
        clear zw
    end
end

% We will solve for pressure, water saturation (oil saturation follows via
% the definition of saturations), polymer concentration and well rates +
% bhp.
primaryVars = {'pressure', 'sW', 'polymer', wellVarNames{:}};

% Evaluate relative permeability
sO  = 1 - sW;
sO0 = 1 - sW0;

[krW, krO] = model.evaluateRelPerm({sW, sO});

% Multipliers for properties
[pvMult, transMult, mobMult, pvMult0] = getMultipliers(model.fluid, p, p0);

% Modifiy relperm by mobility multiplier (if any)
krW = mobMult.*krW; krO = mobMult.*krO;

% Compute transmissibility
T = op.T.*transMult;

% Gravity contribution
gdz = model.getGravityGradient();

% Evaluate water and polymer properties
ads  = effads(c, cmax, model);
ads0 = effads(c0, cmax0, model);
[vW, vP, bW, ~, mobW, mobP, rhoW, pW, upcw, a] = ...
    getFluxAndPropsWaterPolymer_BO(model, p, sW, c, ads, ...
    krW, T, gdz);
bW0 = model.fluid.bW(p0);

% Evaluate oil properties
[vO, bO, mobO, rhoO, p, upco] = getFluxAndPropsOil_BO(model, p, sO, krO, T, gdz);
bO0 = getbO_BO(model, p0);

if model.outputFluxes
    state = model.storeFluxes(state, vW, vO, vP);
end
if model.extraStateOutput
    state = model.storebfactors(state, bW, bO, []);
    state = model.storeMobilities(state, mobW, mobO, mobP);
    state = model.storeUpstreamIndices(state, upcw, upco, []);
end

% EQUATIONS ---------------------------------------------------------------
% Upstream weight b factors and multiply by interface fluxes to obtain the
% fluxes at standard conditions.
bOvO = op.faceUpstr(upco, bO).*vO;
bWvW = op.faceUpstr(upcw, bW).*vW;
bWvP = op.faceUpstr(upcw, bW).*vP;

% Conservation of mass for water
water = (op.pv/dt).*( pvMult.*bW.*sW - pvMult0.*bW0.*sW0 ) + op.Div(bWvW);

% Conservation of mass for oil
oil = (op.pv/dt).*( pvMult.*bO.*sO - pvMult0.*bO0.*sO0 ) + op.Div(bOvO);

% Conservation of polymer in water:
poro = model.rock.poro;
f    = model.fluid;
polymer = (op.pv.*(1-f.dps)/dt).*(pvMult.*bW.*sW.*c - ...
   pvMult0.*bW0.*sW0.*c0) + (op.pv/dt).* ...
   ( f.rhoR.*((1-poro)./poro).*(ads-ads0) ) + op.Div(bWvP);

eqs   = {water, oil, polymer};
names = {'water', 'oil', 'polymer'};
types = {'cell', 'cell', 'cell'};

% Add in any fluxes / source terms prescribed as boundary conditions.
rho = {rhoW, rhoO};
mob = {mobW, mobO};
sat = {sW, sO};
[eqs, qBC, qRes, BCTocellMap, qSRC, srcCells] = addFluxesFromSourcesAndBC(...
   model, eqs, {pW, p}, rho, mob,  ...
   sat, drivingForces);

if model.outputFluxes
    state = model.storeBoundaryFluxes(state, qRes{1}, qRes{2}, [], drivingForces);
end
% Add polymer boundary conditions
if ~isempty(drivingForces.bc) && isfield(drivingForces.bc, 'poly')
   injInx  = qBC{1} > 0; % water inflow indices
   cbc     = (BCTocellMap')*c; % BCTocellMap' = cellToBCMap
   cbc(injInx) = drivingForces.bc.poly(injInx);
   eqs{3}  = eqs{3} - BCTocellMap*(cbc.*qBC{1});
end

% Add polymer source
if ~isempty(drivingForces.src) && isfield(drivingForces.src, 'poly')
   injInx  = qSRC{1} > 0; % water inflow indices
   csrc    = c(srcCells);
   csrc(injInx) = drivingForces.src.poly(injInx);
   eqs{3}(srcCells) = eqs{3}(srcCells) - csrc.*qSRC{1};
end

% Finally, add in and setup well equations
if ~isempty(W) 
    [eqs, names, types, state.wellSol] = model.insertWellEquations(eqs, names, types, wellSol0, wellSol, qWell, pBH, wellVars, wellMap, p, mob, rho, {}, {c}, dt, opt);
end
problem = LinearizedProblem(eqs, types, names, primaryVars, state, dt);
end


%--------------------------------------------------------------------------


% Effective adsorption, depending of desorption or not
function y = effads(c, cmax, model)
   if model.fluid.adsInx == 2
      y = model.fluid.ads(max(c, cmax));
   else
      y = model.fluid.ads(c);
   end
end



