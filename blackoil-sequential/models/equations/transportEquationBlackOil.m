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

s = model.operators;
f = model.fluid;

disgas = model.disgas;
vapoil = model.vapoil;

% Properties at current timestep
[p, sW, sG, rs, rv, wellSol] = model.getProps(state, ...
                                'pressure', 'water', 'gas', 'rs', 'rv', 'wellSol');
% Properties at previous timestep
[p0, sW0, sG0, rs0, rv0] = model.getProps(state0, ...
                                'pressure', 'water', 'gas', 'rs', 'rv');
% If timestep has been split relative to pressure, linearly interpolate in
% pressure.
if isfield(state, 'timestep')
    dt_frac = dt/state.timestep;
    p = p.*dt_frac + p0.*(1-dt_frac);
end
%Initialization of primary variables ----------------------------------
st  = getCellStatusVO(model, state,  1-sW-sG,   sW,  sG);
st0 = getCellStatusVO(model, state0, 1-sW0-sG0, sW0, sG0);
if ~opt.resOnly,
    if ~opt.reverseMode,
        % define primary varible x and initialize
        x = st{1}.*rs + st{2}.*rv + st{3}.*sG;

        [sW, x] = initVariablesADI(sW, x);

        % define sG, rs and rv in terms of x
        sG = st{2}.*(1-sW) + st{3}.*x;

        if disgas
            rsSat = f.rsSat(p);
            rs = (~st{1}).*rsSat + st{1}.*x;
        end
        if vapoil
            rvSat = f.rvSat(p);
            rv = (~st{2}).*rvSat + st{2}.*x;
        end
    else
        assert(0, 'Backwards solver not supported for splitting');
    end
end
if disgas || vapoil
    gvar = 'x';
else
    gvar = 'sG';
end
primaryVars = {'sW', gvar};

% Evaluate relative permeability
sO  = 1 - sW  - sG;
sO0 = 1 - sW0 - sG0;
[krW, krO, krG] = model.evaluateRelPerm({sW, sO, sG});

% Multipliers for properties
[pvMult, transMult, mobMult, pvMult0] = getMultipliers(model.fluid, p, p0);

% Modifiy relperm by mobility multiplier (if any)
krW = mobMult.*krW; krO = mobMult.*krO; krG = mobMult.*krG;

% Compute transmissibility
T = s.T.*transMult;

% Gravity gradient per face
gdz = model.getGravityGradient();

% Evaluate water properties
[vW, bW, mobW, rhoW, pW, upcw, dpW] = getFluxAndPropsWater_BO(model, p, sW, krW, T, gdz);
bW0 = f.bW(p0);

% Evaluate oil properties
[vO, bO, mobO, rhoO, p, upco, dpO] = getFluxAndPropsOil_BO(model, p, sO, krO, T, gdz, rs, ~st{1});
bO0 = getbO_BO(model, p0, rs0, ~st0{1});

% Evaluate gas properties
bG0 = getbG_BO(model, p0, rv0, ~st0{2});
[vG, bG, mobG, rhoG, pG, upcg, dpG] = getFluxAndPropsGas_BO(model, p, sG, krG, T, gdz, rv, ~st{2});

% Get total flux from state
flux = sum(state.flux, 2);
vT = flux(model.operators.internalConn);

% Sat dependent pressure terms
gp = s.Grad(p);
Gw = gp - dpW;
Go = gp - dpO;
Gg = gp - dpG;

flag = getSaturationUpwind(model.upwindType, state, {Gw, Go, Gg}, vT, s.T, {mobW, mobO, mobG}, s.faceUpstr);

upcw  = flag(:, 1);
upco  = flag(:, 2);
upcg  = flag(:, 3);

% Upstream weighted face mobilities
mobWf = s.faceUpstr(upcw, mobW);
mobOf = s.faceUpstr(upco, mobO);
mobGf = s.faceUpstr(upcg, mobG);
% Tot mob
totMob = mobOf + mobWf + mobGf;

f_w = mobWf./totMob;
f_o = mobOf./totMob;
f_g = mobGf./totMob;

vW = f_w.*(vT + s.T.*mobOf.*(Gw - Go) + s.T.*mobGf.*(Gw - Gg));
vO = f_o.*(vT + s.T.*mobWf.*(Go - Gw) + s.T.*mobGf.*(Go - Gg));
vG = f_g.*(vT + s.T.*mobWf.*(Gg - Gw) + s.T.*mobOf.*(Gg - Go));

bWvW = s.faceUpstr(upcw, bW).*vW;
bOvO = s.faceUpstr(upco, bO).*vO;
bGvG = s.faceUpstr(upcg, bG).*vG;

if disgas
    rsbOvO = s.faceUpstr(upco, rs).*bOvO;
end

if vapoil;
    rvbGvG = s.faceUpstr(upcg, rv).*bGvG;
end

if model.extraStateOutput
    state = model.storebfactors(state, bW, bO, bG);
    state = model.storeMobilities(state, mobW, mobO, mobG);
end
% well equations
if ~isempty(W)
    wflux = vertcat(wellSol.flux);
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

    f_w_w = mobWw./totMobw;
    f_o_w = mobOw./totMobw;
    f_g_w = mobGw./totMobw;


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
        wflux_G = wflux_G + bOqO.*rs(wc);
    end

    if vapoil
        wflux_O = wflux_O + bGqG.*rv(wc);
    end

    for i = 1:numel(W)
        perfind = perf2well == i;
        state.wellSol(i).qOs = sum(double(wflux_O(perfind)));
        state.wellSol(i).qWs = sum(double(wflux_W(perfind)));
        state.wellSol(i).qGs = sum(double(wflux_G(perfind)));
    end
end

[eqs, names, types] = deal(cell(1,2));
if opt.solveForWater
    % water eq:
    wat = (s.pv/dt).*( pvMult.*bW.*sW - pvMult0.*bW0.*sW0 ) + s.Div(bWvW);
    if ~isempty(W)
        wat(wc) = wat(wc) - wflux_W;
    end
else
    wat = [];
end

if opt.solveForGas
    % gas eq:
    if disgas
        gas = (s.pv/dt).*( pvMult.* (bG.* sG  + rs.* bO.*sO) - ...
                              pvMult0.*(bG0.*sG0 + rs0.*bO0.*sO0 ) ) + ...
                 s.Div(bGvG + rsbOvO);
    else
        gas = (s.pv/dt).*( pvMult.*bG.*sG - pvMult0.*bG0.*sG0 ) + s.Div(bGvG);
    end
    if ~isempty(W)
        gas(wc) = gas(wc) - wflux_G;
    end
else
    gas = [];
end

if opt.solveForOil
    % oil eq:
    if vapoil
        oil = (s.pv/dt).*( pvMult.* (bO.* sO  + rv.* bG.* sG) - ...
                              pvMult0.*(bO0.*sO0 + rv0.*bG0.*sG0) ) + ...
                 s.Div(bOvO + rvbGvG);
    else
        oil = (s.pv/dt).*( pvMult.*bO.*sO - pvMult0.*bO0.*sO0 ) + s.Div(bOvO);
    end
    if ~isempty(W)
        oil(wc) = oil(wc) - wflux_O;
    end
else
    oil = [];
end

phaseEqs = {wat, oil, gas};
% Add in any fluxes / source terms prescribed as boundary conditions.
phaseEqs = addFluxesFromSourcesAndBC(model, phaseEqs, ...
                                       {pW, p, pG},...
                                       {rhoW, rhoO, rhoG},...
                                       {mobW, mobO, mobG}, ...
                                       {bW, bO, bG},  ...
                                       {sW, sO, sG}, ...
                                       drivingForces);
ix = 1;
active = [opt.solveForWater, opt.solveForOil, opt.solveForGas];
enames = {'water', 'oil', 'gas'};
for i = 1:numel(active)
    if active(i)
        names{ix} = enames{i};
        types{ix} = 'cell';
        eqs{ix} = phaseEqs{i};
        if ~model.useCNVConvergence
            eqs{ix} = eqs{ix}.*(dt./s.pv);
        end
        ix = ix + 1;
    end
end
problem = LinearizedProblem(eqs, types, names, primaryVars, state, dt);
end

