function checkComponentMassBalance(model, state0, states, schedule, n)
% Check mass balance of a simulator run and print to screen
    if nargin < 5
        n = numel(states);
    end
    
    if ~isfield(state0, 'L')
        state0 = model.computeFlash(state0, schedule.step.val(1), 1);
    end
    
    states = states(1:n);
    dt = schedule.step.val(1:n);

    ws = cellfun(@(x) x.wellSol, states, 'UniformOutput', false);    
    ncomp = numel(states{1}.components);
    
    fracx = getMassFraction(states{end}.x, model.EOSModel.fluid);
    fracy = getMassFraction(states{end}.y, model.EOSModel.fluid);

    fracx0 = getMassFraction(state0.x, model.EOSModel.fluid);
    fracy0 = getMassFraction(state0.y, model.EOSModel.fluid);
    
    rhoO0 = model.EOSModel.computeDensity(state0.pressure, state0.x, state0.Z_L, state0.T);
    rhoG0 = model.EOSModel.computeDensity(state0.pressure, state0.y, state0.Z_V, state0.T);
    for i = 1:numel(states)
        if ~isfield(states{i}, 'rho')
            states{i} = computeDensities(model, states{i});
        end
    end
    watOffset = model.water;
    
    [oilMass, oilMass0, gasMass, gasMass0, injMass, prodMass] = deal(0);
    for i = 1:ncomp
        wcomp = bsxfun(@times, getWellComponent(ws, i), dt);
        wcomp = wcomp(:);
        
        injected = sum(wcomp(wcomp > 0));
        produced = abs(sum(wcomp(wcomp < 0)));
        
        oil = sum(fracx(:, i).*states{end}.s(:, 1 + watOffset).*states{end}.rho(:, 1 + watOffset).*model.operators.pv);
        gas = sum(fracy(:, i).*states{end}.s(:, 2 + watOffset).*states{end}.rho(:, 2 + watOffset).*model.operators.pv);
        
        oil0 = sum(fracx0(:, i).*state0.s(:, 1 + watOffset).*rhoO0.*model.operators.pv);
        gas0 = sum(fracy0(:, i).*state0.s(:, 2 + watOffset).*rhoG0.*model.operators.pv);
        
        printTable(model.EOSModel.fluid.names{i}, oil0, oil, gas0, gas, injected, produced)

        oilMass = oilMass + oil;
        gasMass = gasMass + gas;
        oilMass0 = oilMass0 + oil0;
        gasMass0 = gasMass0 + gas0;
        injMass = injMass + injected;
        prodMass = prodMass + produced;
    end
    printTable([], oilMass0, oilMass, gasMass0, gasMass, injMass, prodMass)
end

function printTable(name, oil0, oil, gas0, gas, injected, produced)
    doil = oil - oil0;
    dgas = gas - gas0;
    if isempty(name)
        fprintf('Total mass:\n')
    else
        fprintf('Component ''%s'':\n', name);
    end
    fprintf('* Oil mass: From %1.2e -> %1.2e (net %+1.2e)\n', oil0, oil, doil);
    fprintf('* Gas mass: From %1.2e -> %1.2e (net %+1.2e)\n', gas0, gas, dgas);
    fprintf('* Injected: %1.2e\n', injected);
    fprintf('* Produced: %1.2e\n', produced);
    
    in = injected + oil0 + gas0;
    out = produced + oil + gas;
    fprintf('* Start: %1.2e, End %1.2e - %1.2f%%.\n', in, out, 100*in./out);
end

function frac = getMassFraction(components, fluid)
    ncomp = numel(components);
    nv = numel(components{1});
    mass = zeros(nv, ncomp);
    for i = 1:numel(components)
        mass(:, i) = fluid.molarMass(i).*components{i};
    end
    frac = bsxfun(@rdivide, mass, sum(mass, 2));
end

function d = getWellComponent(ws, compix)
    nw = numel(ws{1});
    d = zeros(numel(ws), nw);
    for i = 1:nw
        d(:, i) = cellfun(@(w) w(i).components(compix), ws);
    end
end

function state = computeDensities(model, state)
    rhoO = model.EOSModel.computeDensity(state.pressure, state.x, state.Z_L, state.T);
    rhoG = model.EOSModel.computeDensity(state.pressure, state.y, state.Z_V, state.T);
    
    if model.water
        state.rho = [model.fluid.rhoWS.*model.fluid.bW(state.pressure), rhoO, rhoG];
    else
        state.rho = [rhoO, rhoG];
    end
end
%{
Copyright 2009-2017 SINTEF ICT, Applied Mathematics.

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
