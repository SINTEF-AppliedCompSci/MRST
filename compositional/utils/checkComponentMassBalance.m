function checkComponentMassBalance(model, state0, states, schedule, n)
% Check mass balance of a simulator run and print to screen

    if nargin < 5
        n = numel(states);
    end
    
    if ~isfield(state0, 'x')
        state0 = model.computeFlash(state0, schedule.step.val(1), 1);
    end
    
    states = states(1:n);
    dt = schedule.step.val(1:n);

    ws = cellfun(@(x) x.wellSol, states, 'UniformOutput', false);    
    ncomp = numel(model.getComponentNames());
    
    if isfield(states{end}, 'xMass')
        fracx = states{end}.xMass;
    else
        fracx = getMassFraction(states{end}.x, model.EOSModel.fluid);
    end
    if isfield(states{end}, 'yMass')
        fracy = states{end}.yMass;
    else
        fracy = getMassFraction(states{end}.y, model.EOSModel.fluid);
    end

    fracx0 = getMassFraction(state0.x, model.EOSModel.fluid);
    fracy0 = getMassFraction(state0.y, model.EOSModel.fluid);
    oilIx = 1 + model.water;
    gasIx = 1 + model.water + model.oil;

    if isfield(state0, 'rho')
        rhoO0 = state0.rho(:, oilIx);
        rhoG0 = state0.rho(:, gasIx);
    else
        rhoO0 = model.EOSModel.PropertyModel.computeDensity(state0.pressure, state0.x, state0.Z_L, state0.T);
        rhoG0 = model.EOSModel.PropertyModel.computeDensity(state0.pressure, state0.y, state0.Z_V, state0.T);
    end
    for i = 1:numel(states)
        if ~isfield(states{i}, 'rho')
            states{i} = computeDensities(model, states{i});
        end
    end
    pv = sum(states{end}.s, 2).*model.operators.pv;
    if isfield(model.fluid', 'pvMultR')
        pv = pv.*model.fluid.pvMultR(states{end}.pressure);
    end
    pv0 = sum(state0.s, 2).*model.operators.pv;
    if isfield(model.fluid', 'pvMultR')
        pv0 = pv0.*model.fluid.pvMultR(state0.pressure);
    end
    [oilMass, oilMass0, gasMass, gasMass0, injMass, prodMass] = deal(0);
    
    for i = 1:ncomp
        wcomp = bsxfun(@times, getWellComponent(ws, i), dt);
        wcomp = wcomp(:);
        
        injected = sum(wcomp(wcomp > 0));
        produced = abs(sum(wcomp(wcomp < 0)));
        

        
        
        oil = sum(fracx(:, i).*states{end}.s(:, oilIx).*states{end}.rho(:, oilIx).*pv);
        gas = sum(fracy(:, i).*states{end}.s(:, gasIx).*states{end}.rho(:, gasIx).*pv);
        
        oil0 = sum(fracx0(:, i).*state0.s(:, oilIx).*rhoO0.*pv0);
        gas0 = sum(fracy0(:, i).*state0.s(:, gasIx).*rhoG0.*pv0);
        
        if nargout == 0
            printTable(model.EOSModel.fluid.names{i}, oil0, oil, gas0, gas, injected, produced)
        end

        oilMass = oilMass + oil;
        gasMass = gasMass + gas;
        oilMass0 = oilMass0 + oil0;
        gasMass0 = gasMass0 + gas0;
        injMass = injMass + injected;
        prodMass = prodMass + produced;
    end
    
    if model.water
        rhoW0 = model.fluid.bW(state0.pressure).*model.fluid.rhoWS;
        rhoW = model.fluid.bW(states{end}.pressure).*model.fluid.rhoWS;
        
        watMass = sum(pv.*states{end}.s(:, 1).*rhoW);
        watMass0 = sum(pv0.*state0.s(:, 1).*rhoW0);
        
        wrate = bsxfun(@times, getWellOutput(ws, 'qWs'), dt);
        wrate = wrate(:).*model.fluid.rhoWS;
        injWat  = sum(wrate(wrate > 0));
        prodWat = abs(sum(wrate(wrate < 0)));
        printTableWater(watMass0, watMass, injWat, prodWat)
    end
    if nargout == 0
        printTable([], oilMass0, oilMass, gasMass0, gasMass, injMass, prodMass)
    end
end

function printTableWater(wat0, wat, injected, produced)
    dwat = wat - wat0;
    fprintf('Water mass:\n')

    fn = @(num) formatMassString(num, 1);
    fprintf('* Water: From %s -> %s (net %s)\n', fn(wat0), fn(wat), fn(dwat));
    fprintf('* Injected: %s\n', fn(injected));
    fprintf('* Produced: %s\n', fn(produced));
    
    in = injected + wat0;
    out = produced + wat;
    fprintf('* Start: %s, End %s - %1.2f%%.\n', fn(in), fn(out), 100*in./out);
end

function printTable(name, oil0, oil, gas0, gas, injected, produced)
    doil = oil - oil0;
    dgas = gas - gas0;
    if isempty(name)
        fprintf('Total mass:\n')
    else
        fprintf('Component ''%s'':\n', name);
    end
    
    fn = @(num) formatMassString(num, 1);
    
    fprintf('* Oil mass: From %s -> %s (net %s)\n', fn(oil0), fn(oil), fn(doil));
    fprintf('* Gas mass: From %s -> %s (net %s)\n', fn(gas0), fn(gas), fn(dgas));
    fprintf('* Injected: %s\n', fn(injected));
    fprintf('* Produced: %s\n', fn(produced));
    
    in = injected + oil0 + gas0;
    out = produced + oil + gas;
    fprintf('* Start: %s, End %s - %1.2f%%.\n', fn(in), fn(out), 100*in./out);
end

function frac = getMassFraction(components, fluid)
    [nv, ncomp] = size(components);
    mass = zeros(nv, ncomp);
    for i = 1:ncomp
        mass(:, i) = fluid.molarMass(i).*components(:, i);
    end
    frac = bsxfun(@rdivide, mass, sum(mass, 2));
end

function d = getWellComponent(ws, compix)
    nw = numel(ws{1});
    d = zeros(numel(ws), nw);
    for i = 1:nw
        d(:, i) = cellfun(@(w) sum(w(i).components(:, compix)), ws);
    end
end

function state = computeDensities(model, state)
    if ~isfield(state, 'Z_L')
        eos = model.EOSModel;
        [Si_L, Si_V, A_L, A_V, B_L, B_V, Bi] = eos.getMixtureFugacityCoefficients(state.pressure, state.T, state.x, state.y, eos.fluid.acentricFactors);
        state.Z_L = eos.computeCompressibilityZ(state.pressure, state.x, A_L, B_L, Si_L, Bi, true);
        state.Z_V = eos.computeCompressibilityZ(state.pressure, state.y, A_V, B_V, Si_V, Bi, false);
    end
    rhoO = model.EOSModel.PropertyModel.computeDensity(state.pressure, state.x, state.Z_L, state.T);
    rhoG = model.EOSModel.PropertyModel.computeDensity(state.pressure, state.y, state.Z_V, state.T);
    
    if model.water
        state.rho = [model.fluid.rhoWS.*model.fluid.bW(state.pressure), rhoO, rhoG];
    else
        state.rho = [rhoO, rhoG];
    end
end
%{
Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.

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
