function [acc, div, names, types] = conservationEquationsBlackOil(state0, state, model, dt, drivingForces)
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
[rs0, rv0, sO0, sG0] = model.getProps(state0, 'rs', 'rv', 'so', 'sg');

if model.water
    sW = model.getProp(state, 'sw');
    sW0 = model.getProp(state0, 'sw');
else
    [sW, sW0] = deal(0);
end

s = model.operators;


[rho, mob, b, pv, p_phase] = model.getProps(state, 'Density', 'Mobility', 'ShrinkageFactors', 'PoreVolume', 'PhasePressures');
[b0, pv0] = model.getProps(state0, 'ShrinkageFactors', 'PoreVolume');


[bW, bO, bG] = deal(b{:});
[bW0, bO0, bG0] = deal(b0{:});

[ff, flags] = model.getProps(state, 'PhaseFlux',  'PhaseUpwindFlag');

nph = 3;
v = cell(nph, 1);
for i = 1:nph
    v{i} = s.faceUpstr(flags{i}, b{i}).*ff{i};
end
[bWvW, bOvO, bGvG] = deal(v{:});

% The first equation is the conservation of the water phase. This equation is
% straightforward, as water is assumed to remain in the aqua phase in the
% black oil model.
if model.water
    water = (1/dt).*( pv.*bW.*sW - pv0.*bW0.*sW0 );
    divWater = s.Div(bWvW);
end

% Second equation: mass conservation equation for the oil phase at surface
% conditions. This is any liquid oil at reservoir conditions, as well as
% any oil dissolved into the gas phase (if the model has vapoil enabled).
if model.vapoil
    % The model allows oil to vaporize into the gas phase. The conservation
    % equation for oil must then include the fraction present in the gas
    % phase.
    rvbGvG = s.faceUpstr(flags{3}, rv).*bGvG;
    % Final equation
    oil = (1/dt).*( pv .* (bO.* sO  + rv.* bG.* sG) - ...
                    pv0.* (bO0.*sO0 + rv0.*bG0.*sG0));
    divOil = s.Div(bOvO + rvbGvG);
else
    oil = (1/dt).*( pv.*bO.*sO - pv0.*bO0.*sO0 );
    divOil = s.Div(bOvO);
end

% Conservation of mass for gas. Again, we have two cases depending on
% whether the model allows us to dissolve the gas phase into the oil phase.
if model.disgas
    % The gas transported in the oil phase.
    rsbOvO = s.faceUpstr(flags{2}, rs).*bOvO;
    
    gas = (1/dt).*( pv.* (bG.* sG  + rs.* bO.* sO) - ...
                    pv0.*(bG0.*sG0 + rs0.*bO0.*sO0 ));
    divGas = s.Div(bGvG + rsbOvO);
else
    gas = (1/dt).*( pv.*bG.*sG - pv0.*bG0.*sG0 );
    divGas = s.Div(bGvG);
end

% Put the set of equations into cell arrays along with their names/types.
if model.water
    acc      = {water, oil, gas};
    div = {divWater, divOil, divGas};
    names    = {'water', 'oil', 'gas'};
    types    = {'cell', 'cell', 'cell'};
else
    acc      = {oil, gas};
    div = {divOil, divGas};
    names    = {'oil', 'gas'};
    types    = {'cell', 'cell'};
end

end

