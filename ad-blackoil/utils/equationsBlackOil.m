function [problem, state] = equationsBlackOil(state0, state, model, dt, drivingForces, varargin)
% Generate linearized problem for a volatile 3Ph system (wet-gas, live-oil).
opt = struct('Verbose',     mrstVerbose,...
    'reverseMode', false,...
    'resOnly',     false,...
    'iteration',   -1);

opt = merge_options(opt, varargin{:});

W = drivingForces.Wells;
assert(isempty(drivingForces.bc) && isempty(drivingForces.src))

% Operators, grid and fluid model.
s = model.operators;
G = model.G;
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

% Pressure dependent transmissibility multiplier 
[trMult, pvMult, pvMult0, transMult] = deal(1);
if isfield(f, 'tranMultR')
    trMult = f.tranMultR(p);
end
% Pressure dependent pore volume multiplier
if isfield(f, 'pvMultR')
    pvMult =  f.pvMultR(p);
    pvMult0 = f.pvMultR(p0);
end
if isfield(f, 'transMult')
   transMult = f.transMult(p); 
end

% Check for capillary pressure (p_cow)
pcOW = 0;
if isfield(f, 'pcOW')
    pcOW  = f.pcOW(sW);
end
%C heck for capillary pressure (p_cog)
pcOG = 0;
if isfield(f, 'pcOG')
    pcOG  = f.pcOG(sG);
end

% Gravity contribution, assert that it is aligned with z-dir
grav = gravity();
grav = grav(1:G.griddim);
%assert(grav(1) == 0 && grav(2) == 0);
%g  = norm(grav);
%dz = s.Grad(G.cells.centroids(:,3));
gdz = s.Grad(G.cells.centroids) * grav';
% Compute transmissibility
T = s.T.*transMult;

% Evaluate relative permeability
sO  = 1 - sW  - sG;
sO0 = 1 - sW0 - sG0;

[krW, krO, krG] = model.evaluteRelPerm({sW, sO, sG});

% Water props (calculated at oil pressure)
bW     = f.bW(p);
bW0 = f.bW(p0);
rhoW   = bW.*f.rhoWS;
% rhoW on face, avarge of neighboring cells (E100, not E300)
rhoWf  = s.faceAvg(rhoW);
mobW   = trMult.*krW./f.muW(p);
dpW    = s.Grad(p-pcOW) - rhoWf.*gdz;
% water upstream-index
upcw  = (double(dpW)<=0);
vW = - s.faceUpstr(upcw, mobW).*T.*dpW;
bWvW = s.faceUpstr(upcw, bW).*vW;
if any(bW < 0)
    warning('Negative water compressibility present!')
end

% Oil props
if disgas
    bO  = f.bO(p,  rs, ~st{1});
    bO0 = f.bO(p0, rs0, ~st0{1}); 
    muO = f.muO(p, rs, ~st{1});
else
    bO  = f.bO(p);
    bO0 = f.bO(p0);
    muO = f.muO(p);
end
if any(bO < 0)
    warning('Negative oil compressibility present!')
end
rhoO   = bO.*(rs*f.rhoGS + f.rhoOS);
rhoOf  = s.faceAvg(rhoO);
mobO   = trMult.*krO./muO;
dpO    = s.Grad(p) - rhoOf.*gdz;
% oil upstream-index
upco = (double(dpO)<=0);
vO   = - s.faceUpstr(upco, mobO).*T.*dpO;
bOvO   = s.faceUpstr(upco, bO).*vO;
if disgas
    rsbOvO = s.faceUpstr(upco, rs).*bOvO;
end

% Gas props (calculated at oil pressure)
if vapoil
    bG  = f.bG(p, rv, ~st{2});
    bG0 = f.bG(p0, rv0, ~st0{2});
    muG = f.muG(p, rv, ~st{2});
else
    bG  = f.bG(p);
    bG0 = f.bG(p0);
    muG = f.muG(p);
end
if any(bG < 0)
    warning('Negative gas compressibility present!')
end
rhoG   = bG.*(rv*f.rhoOS + f.rhoGS);
rhoGf  = s.faceAvg(rhoG);
mobG   = trMult.*krG./muG;
dpG    = s.Grad(p+pcOG) - rhoGf.*gdz;
% gas upstream-index
upcg    = (double(dpG)<=0);
vG = - s.faceUpstr(upcg, mobG).*T.*dpG;
bGvG   = s.faceUpstr(upcg, bG).*vG;
if vapoil
    rvbGvG = s.faceUpstr(upcg, rv).*bGvG;
end

if model.outputFluxes
    state = model.storeFluxes(state, vW, vO, vG);
end

if model.extraStateOutput
    state = model.storebfactors(state, bW, bO, bG);
    state = model.storeMobilities(state, mobW, mobO, mobG);
    state = model.storeUpstreamIndices(state, upcw, upco, upcg);
end

% EQUATIONS -----------------------------------------------------------

% oil eq:
names{1} = 'oil';
if vapoil
    eqs{1} = (s.pv/dt).*( pvMult.* (bO.* sO  + rv.* bG.* sG) - ...
        pvMult0.*(bO0.*sO0 + rv0.*bG0.*sG0) ) + ...
        s.Div(bOvO + rvbGvG);
else
    eqs{1} = (s.pv/dt).*( pvMult.*bO.*sO - pvMult0.*bO0.*sO0 ) + s.Div(bOvO);
end


% water eq:
names{2} = 'water';
eqs{2} = (s.pv/dt).*( pvMult.*bW.*sW - pvMult0.*bW0.*sW0 ) + s.Div(bWvW);

% gas eq:
names{3} = 'gas';
if disgas
    eqs{3} = (s.pv/dt).*( pvMult.* (bG.* sG  + rs.* bO.* sO) - ...
        pvMult0.*(bG0.*sG0 + rs0.*bO0.*sO0 ) ) + ...
        s.Div(bGvG + rsbOvO);
else
    eqs{3} = (s.pv/dt).*( pvMult.*bG.*sG - pvMult0.*bG0.*sG0 ) + s.Div(bGvG);
end
types = {'cell', 'cell', 'cell'};

wm = WellModel();

% well equations
if ~isempty(W)
    wc    = vertcat(W.cells);
    if ~opt.reverseMode
        nperf = numel(wc);
        pw    = p(wc);
        rhows = [f.rhoWS, f.rhoOS, f.rhoGS];
        bw    = {bW(wc), bO(wc), bG(wc)};
        if ~disgas
           % rs supposed to be scalar in this case
            rsw = ones(nperf,1)*rs; rsSatw = ones(nperf,1)*rsSat; %constants
        else
            rsw = rs(wc); rsSatw = rsSat(wc);
        end
        if ~vapoil
            rvw = ones(nperf,1)*rv; rvSatw = ones(nperf,1)*rvSat; %constants
        else
           % rv supposed to be scalar in this case
            rvw = rv(wc); rvSatw = rvSat(wc);
        end
        rw    = {rsw, rvw};
        rSatw = {rsSatw, rvSatw};
        mw    = {mobW(wc), mobO(wc), mobG(wc)};
        s = {sW, 1 - sW - sG, sG};
        
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
        % Force wells to be ADI variables.
        nw = numel(state0.wellSol);
        zw = double2ADI(zeros(nw,1), p0);
        eqs(4:7) = {zw, zw, zw, zw};
        names(4:7) = {'empty', 'empty', 'empty', 'empty'};
        types(4:7) = {'none', 'none', 'none', 'none'};
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
