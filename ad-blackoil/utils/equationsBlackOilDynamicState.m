function [eqs, names, types] = equationsBlackOilDynamicState(state0, state, model, dt, drivingForces)
% Black-oil (new state)

%{
Copyright 2009-2018 SINTEF ICT, Applied Mathematics.

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
% Shorter names for some commonly used parts of the model and forces.

[p, rs, rv, sO, sG] = model.getProps(state, 'pressure', 'rs', 'rv', 'so', 'sg');
[p0, rs0, rv0, sO0, sG0] = model.getProps(state0, 'pressure', 'rs', 'rv', 'so', 'sg');

if model.water
    sW = model.getProp(state, 'sw');
    sW0 = model.getProp(state0, 'sw');
    sat = {sW, sO, sG};
else
    sat = {sO, sG};
    [sW, sW0] = deal(0);
end

s = model.operators;
f = model.fluid;
% Multipliers for properties
[pvMult, transMult, mobMult, pvMult0] = getMultipliers(model.fluid, p, p0);

% Evaluate relative permeability
if model.water
    [krW, krO, krG] = model.evaluateRelPerm(sat);
    % Mobility multiplier for water
    krW = mobMult.*krW;
else
    [krO, krG] = model.evaluateRelPerm(sat);
end

% Modifiy relperm by mobility multiplier (if any)
krO = mobMult.*krO; krG = mobMult.*krG;

% Compute transmissibility
T = s.T.*transMult;

% Gravity gradient per face
gdz = model.getGravityGradient();

[rho, state] = model.evaluatePropertyFunction(state, 'density');
[mu, state] = model.evaluatePropertyFunction(state, 'viscosity');


% Evaluate water properties
if model.water
    [vW, bW, mobW, rhoW, pW, upcw] = getFluxAndPropsWater_BO(model, p, sW, krW, T, gdz);
    bW0 = f.bW(p0);
else
    [vW, bW, mobW, rhoW, pW, upcw] = deal([]);
end
st  = model.getCellStatusVO(state,  sO, sW, sG);
st0  = model.getCellStatusVO(state0,  sO0, sW0, sG0);

% Evaluate oil properties
[vO, bO, mobO, rhoO, p, upco] = getFluxAndPropsOil_BO(model, p, sO, krO, T, gdz, rs, ~st{1});
bO0 = getbO_BO(model, p0, rs0, ~st0{1});

% Evaluate gas properties
bG0 = getbG_BO(model, p0, rv0, ~st0{2});
[vG, bG, mobG, rhoG, pG, upcg] = getFluxAndPropsGas_BO(model, p, sG, krG, T, gdz, rv, ~st{2});

% Store fluxes / properties for debugging / plotting, if requested.
% if model.outputFluxes
%     state = model.storeFluxes(state, vW, vO, vG);
% end
% if model.extraStateOutput
%     state = model.storebfactors(state, bW, bO, bG);
%     state = model.storeDensity(state, rhoW, rhoO, rhoG);
%     state = model.storeMobilities(state, mobW, mobO, mobG);
%     state = model.storeUpstreamIndices(state, upcw, upco, upcg);
% end

% EQUATIONS -----------------------------------------------------------

% Upstream weight b factors and multiply by interface fluxes to obtain the
% fluxes at standard conditions.
bOvO = s.faceUpstr(upco, bO).*vO;
bGvG = s.faceUpstr(upcg, bG).*vG;

% The first equation is the conservation of the water phase. This equation is
% straightforward, as water is assumed to remain in the aqua phase in the
% black oil model.
if model.water
    bWvW = s.faceUpstr(upcw, bW).*vW;
    water = (s.pv/dt).*( pvMult.*bW.*sW - pvMult0.*bW0.*sW0 );
    divWater = s.Div(bWvW);
end

% Second equation: mass conservation equation for the oil phase at surface
% conditions. This is any liquid oil at reservoir conditions, as well as
% any oil dissolved into the gas phase (if the model has vapoil enabled).
if model.vapoil
    % The model allows oil to vaporize into the gas phase. The conservation
    % equation for oil must then include the fraction present in the gas
    % phase.
    rvbGvG = s.faceUpstr(upcg, rv).*bGvG;
    % Final equation
    oil = (s.pv/dt).*( pvMult .* (bO.* sO  + rv.* bG.* sG) - ...
                       pvMult0.*(bO0.*sO0 + rv0.*bG0.*sG0));
    divOil = s.Div(bOvO + rvbGvG);
else
    oil = (s.pv/dt).*( pvMult.*bO.*sO - pvMult0.*bO0.*sO0 );
    divOil = s.Div(bOvO);
end

% Conservation of mass for gas. Again, we have two cases depending on
% whether the model allows us to dissolve the gas phase into the oil phase.
if model.disgas
    % The gas transported in the oil phase.
    rsbOvO = s.faceUpstr(upco, rs).*bOvO;
    
    gas = (s.pv/dt).*( pvMult.* (bG.* sG  + rs.* bO.* sO) - ...
                       pvMult0.*(bG0.*sG0 + rs0.*bO0.*sO0 ));
    divGas = s.Div(bGvG + rsbOvO);
else
    gas = (s.pv/dt).*( pvMult.*bG.*sG - pvMult0.*bG0.*sG0 );
    divGas = s.Div(bGvG);
end

% Put the set of equations into cell arrays along with their names/types.
if model.water
    eqs      = {water, oil, gas};
    divTerms = {divWater, divOil, divGas};
    names    = {'water', 'oil', 'gas'};
    types    = {'cell', 'cell', 'cell'};
%     rho      = {rhoW, rhoO, rhoG};
%     mob      = {mobW, mobO, mobG};
%     pressures = {pW, p, pG};
else
    eqs      = {oil, gas};
    divTerms = {divOil, divGas};
    names    = {'oil', 'gas'};
    types    = {'cell', 'cell'};
%     rho      = {rhoO, rhoG};
%     mob      = {mobO, mobG};
%     pressures = {p, pG};
end
% Add in any fluxes / source terms prescribed as boundary conditions.
% dissolved = model.getDissolutionMatrix(rs, rv);
% 
% [eqs, state] = addBoundaryConditionsAndSources(model, eqs, names, types, state, ...
%                                                  pressures, sat, mob, rho, ...
%                                                  dissolved, {}, ...
%                                                  drivingForces);


% Finally, adding divergence terms to equations
for i = 1:numel(divTerms)
    eqs{i} = eqs{i} + divTerms{i};
end

end

