function [problem, state] = transportEquationOilWater(state0, state, model, dt, drivingForces, varargin)

opt = struct('Verbose', mrstVerbose, ...
             'reverseMode', false,...
             'scaling', [],...
             'resOnly', false,...
             'history', [],...
             'solveForWater', false, ...
             'solveForOil', true, ...
             'iteration', -1, ...
             'stepOptions', []);  % Compatibility only

opt = merge_options(opt, varargin{:});
W = drivingForces.W;

op = model.operators;
solveAllPhases = opt.solveForWater && opt.solveForOil;
[p, sW, sO, wellSol] = model.getProps(state, 'pressure', 'water', 'oil', 'wellsol');
[sW0, sO0] = model.getProps(state0, 'water', 'oil');

% Initialization of independent variables ----------------------------------
assert(~opt.reverseMode, 'Backwards solver not supported for splitting');
if solveAllPhases
    if ~opt.resOnly
        [sW, sO] = model.AutoDiffBackend.initVariablesAD(sW, sO);
    end
    primaryVars = {'sW', 'sO'};
    sT = sO + sW;
    s = {sW./sT, sO./sT};
    s0 = state.s;
else
    if ~opt.resOnly
        sW = model.AutoDiffBackend.initVariablesAD(sW);
    end
    primaryVars = {'sW'};
    sO = 1 - sW;
    sT = ones(size(value(sW)));
    s = {sW, sO};
end
state = model.setProp(state, 's', s);
state = model.initStateFunctionContainers(state);
% -------------------------------------------------------------------------
[b, pv] = model.getProps(state, 'ShrinkageFactors', 'PoreVolume');
[b0, pv0] = model.getProps(state0, 'ShrinkageFactors', 'PoreVolume');

[mob, pc, rhogdz, rho, T] = model.getProps(state, 'Mobility', ...
                                               'CapillaryPressure', ...
                                               'GravityPotentialDifference', ...
                                               'Density', ...
                                               'Transmissibility');
if solveAllPhases
    state.s = s0;
end
[bW, bO] = deal(b{:});
[bW0, bO0] = deal(b0{:});
[mobW, mobO] = deal(mob{:});
[rhoW, rhoO] = deal(rho{:});

Gw = -rhogdz{1};
Go = -rhogdz{2};
if ~isempty(pc{1})
    Gw = Gw - op.Grad(pc{1});
end

if model.extraStateOutput
    state = model.storebfactors(state, bW, bO, []);
    state = model.storeMobilities(state, mobW, mobO, []);
    state = model.storeDensity(state, rhoW, rhoO, []);
end

if ~isempty(W)
    wflux = sum(vertcat(wellSol.flux), 2);
    perf2well = getPerforationToWellMapping(W);
    wc = vertcat(W.cells);
    
    mobWw = mobW(wc);
    mobOw = mobO(wc);
    totMobw = mobWw + mobOw;

    f_w_w = sT(wc).*mobWw./totMobw;
    f_o_w = sT(wc).*mobOw./totMobw;

    isInj = wflux > 0;
    compWell = vertcat(W.compi);
    compPerf = compWell(perf2well, :);

    f_w_w(isInj) = compPerf(isInj, 1);
    f_o_w(isInj) = compPerf(isInj, 2);

    bWqW = bW(wc).*f_w_w.*wflux;
    bOqO = bO(wc).*f_o_w.*wflux;

    % Store well fluxes
    wflux_O = value(bOqO);
    wflux_W = value(bWqW);
    
    for i = 1:numel(W)
        perfind = perf2well == i;
        state.wellSol(i).qOs = sum(wflux_O(perfind));
        state.wellSol(i).qWs = sum(wflux_W(perfind));
    end

end
% Get total flux from state
flux = sum(state.flux, 2);
vT = flux(model.operators.internalConn);

sat = {sW, sO};
G = {Gw, Go};

components = {{sT.*bW,  []}, ...
              {[],      sT.*bO}};
upstr = model.operators.faceUpstr;
[q_phase, q_components] = computeSequentialFluxes(...
    state, G, vT, T, mob, rho, components, upstr, model.upwindType);
[waterFlux, oilFlux] = deal(q_components{:});

if model.outputFluxes
    state = model.storeFluxes(state, q_phase{:}, []);
end

if opt.solveForWater
    wat = (1/dt).*(pv.*bW.*sW - pv0.*bW0.*sW0);
    if ~isempty(W)
        wat(wc) = wat(wc) - bWqW;
    end
end

if opt.solveForOil
    oil = (1/dt).*(pv.*bO.*sO - pv0.*bO0.*sO0);
    if ~isempty(W)
        oil(wc) = oil(wc) - bOqO;
    end
end

if solveAllPhases
    eqs = {wat, oil};
    fluxes = {waterFlux, oilFlux};
    names = {'water', 'oil'};    
    types = {'cell', 'cell'};
elseif opt.solveForOil
    eqs = {oil};
    fluxes = {oilFlux};
    names = {'oil'};
    types = {'cell'};
else
    eqs = {wat};
    fluxes = {waterFlux};
    names = {'water'};
    types = {'cell'};
end


[eqs, state] = addBoundaryConditionsAndSources(model, eqs, names, types, state, ...
                                     {p, p}, sat, mob, rho, ...
                                     {}, {}, ...
                                     drivingForces);
for i = 1:numel(eqs)
    eqs{i} = op.AccDiv(eqs{i}, fluxes{i});
    if ~model.useCNVConvergence
        eqs{i} = eqs{i}.*(dt./op.pv);
    end
end

problem = LinearizedProblem(eqs, types, names, primaryVars, state, dt);
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
