function [problem, state] = equationsBlackOil(state0, state, model, dt, drivingForces, varargin)
% Generate linearized problem for a volatile 3Ph system (wet-gas, live-oil).
opt = struct('Verbose',     mrstVerbose,...
    'reverseMode', false,...
    'resOnly',     false,...
    'iteration',   -1);

opt = merge_options(opt, varargin{:});

W = drivingForces.Wells;

% Operators, grid and fluid model.
s = model.operators;
f = model.fluid;

% Can gas dissolve into oil phase (Rs)?
disgas = model.disgas;
% Can oil be present as vapor in gas phase (Rv)?
vapoil = model.vapoil;

% Properties at current timestep
[p, sW, sG, rs, rv, wellSol] = model.getProps(state, ...
    'pressure', 'water', 'gas', 'rs', 'rv', 'wellSol');
% Properties at previous timestep
[p0, sW0, sG0, rs0, rv0] = model.getProps(state0, ...
    'pressure', 'water', 'gas', 'rs', 'rv');

bhp    = vertcat(wellSol.bhp);
qWs    = vertcat(wellSol.qWs);
qOs    = vertcat(wellSol.qOs);
qGs    = vertcat(wellSol.qGs);

%Initialization of primary variables ----------------------------------
st  = getCellStatusVO(state,  1-sW-sG,   sW,  sG,  disgas, vapoil);
st0 = getCellStatusVO(state0, 1-sW0-sG0, sW0, sG0, disgas, vapoil);
if ~opt.resOnly,
    if ~opt.reverseMode,
        % define primary varible x and initialize
        x = st{1}.*rs + st{2}.*rv + st{3}.*sG;
        
        [p, sW, x, qWs, qOs, qGs, bhp] = ...
            initVariablesADI(p, sW, x, qWs, qOs, qGs, bhp);
        [sG, rs, rv, rsSat, rvSat] = calculateHydrocarbonsFromStatus(f, st, 1-sW, x, rs, rv, p, disgas, vapoil);
        
    else
        x0 = st0{1}.*rs0 + st0{2}.*rv0 + st0{3}.*sG0;
        
        [p0, sW0, x0, zw, zw, zw, zw] = ...
            initVariablesADI(p0, sW0, x0, ...
            zeros(size(qWs)) , zeros(size(qOs)) , ...
            zeros(size(qGs)) , zeros(size(bhp)));                %#ok
        clear zw
        [sG0, rs0, rv0] = calculateHydrocarbonsFromStatus(f, st0, 1-sW, x0, rs0, rv0, p0, disgas, vapoil);
    end
else
    [sG, rs, rv, rsSat, rvSat] = calculateHydrocarbonsFromStatus(f, st, 1-sW, 0, rs, rv, p, disgas, vapoil);
end
if disgas || vapoil
    % X is either Rs, Rv or Sg, depending on each cell's saturation status
    gvar = 'x';
else
    gvar = 'sG';
end
% We will solve for pressure, water and gas saturation (oil saturation
% follows via the definition of saturations) and well rates + bhp.
primaryVars = {'pressure', 'sW', gvar, 'qWs', 'qOs', 'qGs', 'bhp'};

% Multipliers for properties
[pvMult, transMult, mobMult, pvMult0] = getMultipliers(model.fluid, p, p0);

% Evaluate relative permeability
sO  = 1 - sW  - sG;
sO0 = 1 - sW0 - sG0;
[krW, krO, krG] = model.evaluteRelPerm({sW, sO, sG});
krW = mobMult.*krW;
krO = mobMult.*krO;
krG = mobMult.*krG;

% Gravity contribution
gdz = model.getGravityGradient();
% Compute transmissibility
T = s.T.*transMult;

% Water props
[vW, bW, mobW, rhoW, pW, upcw] = getFluxAndPropsWater_BO(model, p, sW, krW, T, gdz);
bW0 = f.bW(p0);

% Oil props
bO0 = getbO_BO(model, p0, rs0, ~st0{1});
[vO, bO, mobO, rhoO, p, upco] = getFluxAndPropsOil_BO(model, p, sO, krO, T, gdz, rs, ~st{1});

% Gas props
bG0 = getbG_BO(model, p0, rv0, ~st0{2});
[vG, bG, mobG, rhoG, pG, upcg] = getFluxAndPropsGas_BO(model, p, sG, krG, T, gdz, rv, ~st{2});

% Upstream weight b factors and multiply by interface fluxes to obtain the
% fluxes at standard conditions.

bOvO = s.faceUpstr(upco, bO).*vO;
bWvW = s.faceUpstr(upcw, bW).*vW;
bGvG = s.faceUpstr(upcg, bG).*vG;

if model.outputFluxes
    state = model.storeFluxes(state, vW, vO, vG);
end

if model.extraStateOutput
    state = model.storebfactors(state, bW, bO, bG);
    state = model.storeMobilities(state, mobW, mobO, mobG);
    state = model.storeUpstreamIndices(state, upcw, upco, upcg);
end

% EQUATIONS -----------------------------------------------------------

% Conservation of mass for oil
names{1} = 'oil';
if vapoil
    % The model allows oil to vaporize into the gas phase. The conservation
    % equation for oil must then include the fraction present in the gas
    % phase.
    rvbGvG = s.faceUpstr(upcg, rv).*bGvG;
    
    eqs{1} = (s.pv/dt).*( pvMult.* (bO.* sO  + rv.* bG.* sG) - ...
        pvMult0.*(bO0.*sO0 + rv0.*bG0.*sG0) ) + ...
        s.Div(bOvO + rvbGvG);
else
    eqs{1} = (s.pv/dt).*( pvMult.*bO.*sO - pvMult0.*bO0.*sO0 ) + s.Div(bOvO);
end


% Conservation of mass for water
names{2} = 'water';
eqs{2} = (s.pv/dt).*( pvMult.*bW.*sW - pvMult0.*bW0.*sW0 ) + s.Div(bWvW);

% Conservation of mass for gas
names{3} = 'gas';
if disgas
    rsbOvO = s.faceUpstr(upco, rs).*bOvO;

    eqs{3} = (s.pv/dt).*( pvMult.* (bG.* sG  + rs.* bO.* sO) - ...
        pvMult0.*(bG0.*sG0 + rs0.*bO0.*sO0 ) ) + ...
        s.Div(bGvG + rsbOvO);
else
    eqs{3} = (s.pv/dt).*( pvMult.*bG.*sG - pvMult0.*bG0.*sG0 ) + s.Div(bGvG);
end

eqs([2, 1, 3]) = addFluxesFromSourcesAndBC(model, eqs([2, 1, 3]), ...
                                               {pW, p, pG},...
                                               {rhoW,     rhoO, rhoG},...
                                               {mobW,     mobO, mobG}, ...
                                               {bW, bO, bG},  ...
                                               {sW, sO, sG}, ...
                                               drivingForces);

types = {'cell', 'cell', 'cell'};

wm = WellModel();

% well equations
if ~isempty(W)
    wc    = vertcat(W.cells);
    if ~opt.reverseMode
        
        pw    = p(wc);
        rhows = [f.rhoWS, f.rhoOS, f.rhoGS];
        bw    = {bW(wc), bO(wc), bG(wc)};
        
        [rw, rSatw] = wm.getResSatWell(model, wc, rs, rv, rsSat, rvSat);
        mw    = {mobW(wc), mobO(wc), mobG(wc)};
        s = {sW(wc), sO(wc), sG(wc)};
        
        [cqs, weqs, ctrleqs, wc, state.wellSol]  = wm.computeWellFlux(model, W, wellSol, ...
            bhp, {qWs, qOs, qGs}, pw, rhows, bw, mw, s, rw,...
            'maxComponents', rSatw, ...
            'nonlinearIteration', opt.iteration);
        
        eqs(4:6) = weqs;
        eqs{7} = ctrleqs;
        
        eqs{1}(wc) = eqs{1}(wc) - cqs{2}; % Add src to oil eq
        eqs{2}(wc) = eqs{2}(wc) - cqs{1}; % Add src to water eq
        eqs{3}(wc) = eqs{3}(wc) - cqs{3}; % Add src to gas eq
        
        names(4:7) = {'waterWells', 'oilWells', 'gasWells', 'closureWells'};
        types(4:7) = {'perf', 'perf', 'perf', 'well'};
    else
        [eqs(4:7), names(4:7), types(4:7)] = wm.createReverseModeWellEquations(model, state0.wellSol, p0);
    end
end
problem = LinearizedProblem(eqs, types, names, primaryVars, state, dt);
end

function [sG, rs, rv, rsSat, rvSat] = calculateHydrocarbonsFromStatus(fluid, status, sO, x, rs, rv, pressure, disgas, vapoil)
    % define sG, rs and rv in terms of x
    sG = status{2}.*sO + status{3}.*x;
    if disgas
        rsSat = fluid.rsSat(pressure);
        rs = (~status{1}).*rsSat + status{1}.*x;
    else % otherwise rs = rsSat = const
        rsSat = rs;
    end
    if vapoil
        rvSat = fluid.rvSat(pressure);
        rv = (~status{2}).*rvSat + status{2}.*x;
    else % otherwise rv = rvSat = const
        rvSat = rv;
    end
end
