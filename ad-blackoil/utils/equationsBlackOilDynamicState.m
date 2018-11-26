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
else
    [sW, sW0] = deal(0);
end

s = model.operators;
f = model.fluid;
% Multipliers for properties
if 0
    [pvMult, transMult, mobMult, pvMult0] = getMultipliers(model.fluid, p, p0);
else
    [pvMult, transMult, mobMult, pvMult0] = deal(1);
end


% Compute transmissibility
T = s.T.*transMult;

% Gravity gradient per face
gdz = model.getGravityGradient();

props = model.FlowPropertyFunctions;
rho = props.getProperty(model, state, 'Density');
mob = props.getProperty(model, state, 'Mobility');
b = props.getProperty(model, state, 'ShrinkageFactors');
b0 = props.getProperty(model, state0, 'ShrinkageFactors');

[bW, bO, bG] = deal(b{:});
[bW0, bO0, bG0] = deal(b0{:});


p_phase = getPhasePressures(p, state.FlowProps.CapillaryPressure);

nph = 3;
v = cell(nph, 1);
upc = false(numel(double(T)), nph);
for i = 1:nph
%     if i == 2 && model.disgas
%         rhof  = s.faceAvg(rho{i} + rs.*f.rhoGS./f.rhoOS);
%     else
        rhof  = s.faceAvg(rho{i});
%     end
    pot = s.Grad(p_phase{i}) - rhof.*gdz;
    
    upc(:, i) = double(pot)<=0;
    
    dflux = -T.*pot;
    v{i} = dflux.*s.faceUpstr(upc(:, i), b{i}.*(mob{i}));
end

% Upstream weight b factors and multiply by interface fluxes to obtain the
% fluxes at standard conditions.
bWvW = v{1};
bOvO = v{2};
bGvG = v{3};

% The first equation is the conservation of the water phase. This equation is
% straightforward, as water is assumed to remain in the aqua phase in the
% black oil model.
if model.water
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
    rvbGvG = s.faceUpstr(upc(:, 3), rv).*bGvG;
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
    rsbOvO = s.faceUpstr(upc(:, 2), rs).*bOvO;
    
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

