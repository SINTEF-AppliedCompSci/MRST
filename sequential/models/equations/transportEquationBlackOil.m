function [problem, state] = transportEquationBlackOil(state0, state, model, dt, drivingForces, varargin)
% Generate linearized problem for a volatile 3Ph system (wet-gas, live-oil).
opt = struct('Verbose',     mrstVerbose,...
             'reverseMode', false,...
             'resOnly',     false,...
             'solveForWater', false, ...
             'solveForOil', true, ...
             'solveForGas', true, ...
             'iteration',   -1);

opt = merge_options(opt, varargin{:});

W = drivingForces.W;
assert(isempty(drivingForces.bc) && isempty(drivingForces.src))

op = model.operators;
f = model.fluid;

disgas = model.disgas;
vapoil = model.vapoil;

% Properties at current timestep
[p, sW, sG, sO, rs, rv, wellSol] = model.getProps(state, ...
                                'pressure', 'water', 'gas', 'oil', 'rs', 'rv', 'wellSol');
% Properties at previous timestep
[p0, sW0, sG0, sO0, rs0, rv0] = model.getProps(state0, ...
                                'pressure', 'water', 'gas', 'oil', 'rs', 'rv');
% If timestep has been split relative to pressure, linearly interpolate in
% pressure.
if isfield(state, 'timestep')
    dt_frac = dt/state.timestep;
    p = p.*dt_frac + p0.*(1-dt_frac);
end
solveAllPhases = opt.solveForWater && opt.solveForOil && opt.solveForGas;

assert(~opt.reverseMode, 'Backwards solver not supported for splitting');

%Initialization of primary variables ----------------------------------
st  = model.getCellStatusVO(state,  sO,  sW,  sG);
if disgas || vapoil
    x = st{1}.*rs + st{2}.*rv + st{3}.*sG;
    gvar = 'x';
else
    gvar = 'sG';
end

if solveAllPhases
    if ~opt.resOnly
        if disgas || vapoil
            [sW, sO, x] = model.AutoDiffBackend.initVariablesAD(sW, sO, x);
        else
            [sW, sO, sG] = model.AutoDiffBackend.initVariablesAD(sW, sO, sG);
        end
    end
    primaryVars = {'sW', 'sO', gvar};
else
    if ~opt.resOnly
        if disgas || vapoil
            [sW, x] = model.AutoDiffBackend.initVariablesAD(sW, x);
        else
            [sW, sG] = model.AutoDiffBackend.initVariablesAD(sW, sG);
        end
    end
    primaryVars = {'sW', gvar};
    % Evaluate relative permeability
    sO  = 1 - sW  - sG;
    sO0 = 1 - sW0 - sG0;
end
if disgas || vapoil
    % define sG, rs and rv in terms of x
    sG = st{2}.*(1-sW) + st{3}.*x;
    if disgas
        rsSat = f.rsSat(p);
        rs = (~st{1}).*rsSat + st{1}.*x;
        state = model.setProp(state, 'rs', rs);
    end
    if vapoil
        rvSat = f.rvSat(p);
        rv = (~st{2}).*rvSat + st{2}.*x;
        state = model.setProp(state, 'rv', rv);
    end
end

if solveAllPhases
    sT = sO + sW + sG;
    s = {sW./sT, sO./sT, sG./sT};
    s0 = state.s;
else
    sT = ones(size(value(sW)));
    s = {sW, sO, sG};
end
state = model.initStateFunctionContainers(state);
state = model.setProp(state, 's', s);

if ~opt.resOnly
    if ~opt.reverseMode
    else
        assert(0, 'Backwards solver not supported for splitting');
    end
end
[b, pv] = model.getProps(state, 'ShrinkageFactors', 'PoreVolume');
[b0, pv0] = model.getProps(state0, 'ShrinkageFactors', 'PoreVolume');

[mob, pc, rhogdz, rho, T, pPhase] = model.getProps(state, 'Mobility', ...
                                               'CapillaryPressure', ...
                                               'GravityPotentialDifference', ...
                                               'Density', ...
                                               'Transmissibility', ...
                                               'PhasePressures');
if solveAllPhases
    state.s = s0;
end
[bW, bO, bG] = deal(b{:});
[bW0, bO0, bG0] = deal(b0{:});
[mobW, mobO, mobG] = deal(mob{:});
nph = numel(mob);
G = cell(1, nph);
for i = 1:nph
    G{i} = -rhogdz{i};
    if ~isempty(pc{i})
        G{i} = G{i} - op.Grad(pc{i});
    end
end
% Get total flux from state
flux = sum(state.flux, 2);
vT = flux(model.operators.internalConn);


components = {{sT.*bW,  [],         []}, ...
              {[],      sT.*bO,     sT.*rv.*bG}, ...
              {[],      sT.*bO.*rs, sT.*bG}};
upstr = model.operators.faceUpstr;
[q_phase, q_components] = computeSequentialFluxes(...
    state, G, vT, T, mob, rho, components, upstr, model.upwindType);
[waterFlux, oilFlux, gasFlux] = deal(q_components{:});


if model.outputFluxes
    state = model.storeFluxes(state, q_phase{:});
end

if model.extraStateOutput
    state = model.storebfactors(state, bW, bO, bG);
    state = model.storeMobilities(state, mobW, mobO, mobG);
    state = model.storeDensity(state, rhoW, rhoO, rhoG);
end
% well equations
if ~isempty(W)
    wflux = sum(vertcat(wellSol.flux), 2);
    cqs = vertcat(wellSol.cqs);

    perf2well = getPerforationToWellMapping(W);
    wc    = vertcat(W.cells);

    mobWw = mobW(wc);
    mobOw = mobO(wc);
    mobGw = mobG(wc);
    totMobw = mobWw + mobOw + mobGw;


    isInj = wflux > 0;
    compWell = vertcat(W.compi);
    compPerf = compWell(perf2well, :);

    f_w_w = sT(wc).*mobWw./totMobw;
    f_o_w = sT(wc).*mobOw./totMobw;
    f_g_w = sT(wc).*mobGw./totMobw;


    f_w_w(isInj) = compPerf(isInj, 1);
    f_o_w(isInj) = compPerf(isInj, 2);
    f_g_w(isInj) = compPerf(isInj, 3);


    bWqW = bW(wc).*f_w_w.*wflux;
    bOqO = bO(wc).*f_o_w.*wflux;
    bGqG = bG(wc).*f_g_w.*wflux;

    if 1
        bWqW(isInj) = cqs(isInj, 1);
        bOqO(isInj) = cqs(isInj, 2);
        bGqG(isInj) = cqs(isInj, 3);
    end
    % Store well fluxes
    wflux_O = bOqO;
    wflux_W = bWqW;
    wflux_G = bGqG;

    if disgas
        wflux_G = wflux_G + ~isInj.*bOqO.*rs(wc);
    end

    if vapoil
        wflux_O = wflux_O + ~isInj.*bGqG.*rv(wc);
    end

    for i = 1:numel(W)
        perfind = perf2well == i;
        state.wellSol(i).qOs = sum(value(wflux_O(perfind)));
        state.wellSol(i).qWs = sum(value(wflux_W(perfind)));
        state.wellSol(i).qGs = sum(value(wflux_G(perfind)));
    end
end

[eqs, names, types, fluxes] = deal(cell(1,2 + solveAllPhases));
[types{:}] = deal('cell');
ix = 1;
if opt.solveForWater
    % water eq:
    wat = (1/dt).*( pv.*bW.*sW - pv0.*bW0.*sW0 );
    if ~isempty(W)
        wat(wc) = wat(wc) - wflux_W;
    end
    fluxes{ix} = waterFlux;
    eqs{ix} = wat;
    names{ix} = 'water';
    ix = ix + 1;
end

if opt.solveForOil
    % oil eq:
    if vapoil
        oil = (1/dt).*( pv.* (bO.* sO  + rv.* bG.* sG) - ...
                        pv0.*(bO0.*sO0 + rv0.*bG0.*sG0) ) + ...
                 s.Div(oilFlux);
    else
        oil = (1/dt).*( pv.*bO.*sO - pv0.*bO0.*sO0 );
    end
    if ~isempty(W)
        oil(wc) = oil(wc) - wflux_O;
    end
    fluxes{ix} = oilFlux;
    eqs{ix} = oil;
    names{ix} = 'oil';
    ix = ix + 1;
end

if opt.solveForGas
    % gas eq:
    if disgas
        gas = (1/dt).*( pv.* (bG.* sG  + rs.* bO.*sO) - ...
                        pv0.*(bG0.*sG0 + rs0.*bO0.*sO0 ) );
    else
        gas = (1/dt).*( pv.*bG.*sG - pv0.*bG0.*sG0 );
    end
    if ~isempty(W)
        gas(wc) = gas(wc) - wflux_G;
    end
    fluxes{ix} = gasFlux;
    eqs{ix} = gas;
    names{ix} = 'gas';
end
mob = {mobW, mobO, mobG};
sat = {sW, sO, sG};
dissolved = model.getDissolutionMatrix(rs, rv);

[eqs, state] = addBoundaryConditionsAndSources(model, eqs, names, types, state, ...
                                     pPhase, sat, mob, rho, ...
                                     dissolved, {}, ...
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

