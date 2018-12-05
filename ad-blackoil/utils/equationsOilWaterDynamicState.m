function [eqs, names, types] = equationsBlackOilDynamicState(state0, state, model, dt, drivingForces)
% Black-oil (new state)
timer = tic();
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

[p, sW, sO] = model.getProps(state, 'pressure', 'sw', 'so');
[p0, sW0, sO0] = model.getProps(state0, 'pressure', 'sw', 'so');



s = model.operators;

T = s.T;
% Gravity gradient per face
gdz = model.getGravityGradient();

props = model.FlowPropertyFunctions;
rho = props.getProperty(model, state, 'Density');
mob = props.getProperty(model, state, 'Mobility');
b = props.getProperty(model, state, 'ShrinkageFactors');
pv = props.getProperty(model, state, 'PoreVolume');

b0 = props.getProperty(model, state0, 'ShrinkageFactors');
pv0 = props.getProperty(model, state0, 'PoreVolume');

[bW, bO] = deal(b{:});
[bW0, bO0] = deal(b0{:});
p_phase = getPhasePressures(p, state.FlowProps.CapillaryPressure);
nph = 2;
v = cell(nph, 1);
upc = false(numel(double(T)), nph);
for i = 1:nph
    rhof  = s.faceAvg(rho{i});
    pot = s.Grad(p_phase{i}) - rhof.*gdz;
    
    upc(:, i) = double(pot)<=0;
    
    dflux = -T.*pot;
    v{i} = dflux.*s.faceUpstr(upc(:, i), b{i}.*(mob{i}));
end
[bWvW, bOvO] = deal(v{:});

water = (1/dt).*( pv.*bW.*sW - pv0.*bW0.*sW0 );
divWater = s.Div(bWvW);

oil = (1/dt).*( pv.*bO.*sO - pv0.*bO0.*sO0 );
divOil = s.Div(bOvO);
eqs      = {water, oil};
divTerms = {divWater, divOil};
names    = {'water', 'oil'};
types    = {'cell', 'cell'};

for i = 1:numel(divTerms)
    eqs{i} = eqs{i} + divTerms{i};
end

toc(timer);
end

