function [problem, state] = pressureEquationBlackOil(state0, state, model, dt, drivingForces, varargin)

opt = struct('Verbose', mrstVerbose, ...
             'reverseMode', false,...
             'resOnly', false,...
             'redistributeRS', false, ...
             'propsPressure', [], ...
             'iteration', -1);

opt = merge_options(opt, varargin{:});

W = drivingForces.Wells;
assert(isempty(drivingForces.bc) && isempty(drivingForces.src))

s = model.operators;
f = model.fluid;
G = model.G;

disgas = model.disgas;
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


%Initialization of independent variables ----------------------------------
st  = getCellStatusVO(state,  1-sW-sG,   sW,  sG,  disgas, vapoil);
st0 = getCellStatusVO(state0, 1-sW0-sG0, sW0, sG0, disgas, vapoil);
p_prop = opt.propsPressure;
if ~opt.resOnly,
    if ~opt.reverseMode,
        % define primary varible x and initialize
        x = st{1}.*rs + st{2}.*rv + st{3}.*sG;

        [p, qWs, qOs, qGs, bhp] = ...
            initVariablesADI(p, qWs, qOs, qGs, bhp);
        if isempty(p_prop)
            p_prop = p;
        end
        % define sG, rs and rv in terms of x
        sG = st{2}.*(1-sW) + st{3}.*x;
        if disgas
            rsSat = f.rsSat(p_prop);
            rs = (~st{1}).*rsSat + st{1}.*x;
        else % otherwise rs = rsSat = const
            rsSat = rs;
        end
        if vapoil
            rvSat = f.rvSat(p_prop);
            rv = (~st{2}).*rvSat + st{2}.*x;
        else % otherwise rv = rvSat = const
            rvSat = rv;
        end
    else
        assert(0, 'Backwards solver not supported for splitting');
    end
else % resOnly-case compute rsSat and rvSat for use in well eqs
    if isempty(p_prop)
        p_prop = p;
    end
    if disgas, rsSat = f.rsSat(p_prop); else rsSat = rs; end
    if vapoil, rvSat = f.rvSat(p_prop); else rvSat = rv; end
end
sO  = 1- sW  - sG;
sO0 = 1- sW0 - sG0;

if disgas && opt.redistributeRS
    [sG, rs] = redistributeRS(f, p_prop, rs, sG, sO);
end
primaryVars = {'pressure', 'qWs', 'qOs', 'qGs', 'bhp'};


%----------------------------------------------------------------------
%check for p-dependent tran mult:
trMult = 1;
if isfield(f, 'tranMultR'), trMult = f.tranMultR(p_prop); end

%check for p-dependent porv mult:
pvMult = 1; pvMult0 = 1;
if isfield(f, 'pvMultR')
    pvMult =  f.pvMultR(p_prop);
    pvMult0 = f.pvMultR(p0);
end

%check for capillary pressure (p_cow)
pcOW = 0;
if isfield(f, 'pcOW')
    pcOW  = f.pcOW(sW);
end
%check for capillary pressure (p_cog)
pcOG = 0;
if isfield(f, 'pcOG')
    pcOG  = f.pcOG(sG);
end

% FLIUD PROPERTIES ---------------------------------------------------
[krW, krO, krG] = model.evaluteRelPerm({sW, sO, sG});
% Gravity contribution
gdz = model.getGravityGradient();


% WATER PROPS (calculated at oil pressure)
bW     = f.bW(p_prop);
rhoW   = bW.*f.rhoWS;
% rhoW on face, avarge of neighboring cells (E100, not E300)
rhoWf  = s.faceAvg(rhoW);
muW = f.muW(p_prop);
mobW   = trMult.*krW./muW;
dpW    = s.Grad(p-pcOW) - rhoWf.*gdz;

% water upstream-index
upcw  = (double(dpW)<=0);
vW   = -s.faceUpstr(upcw, mobW).*s.T.*dpW;
bWvW =  s.faceUpstr(upcw, bW).*vW;
if any(bW < 0)
    warning('Negative water compressibility present!')
end

% OIL PROPS
if disgas
    bO  = f.bO(p_prop, rs, ~st{1});
    muO = f.muO(p_prop, rs, ~st{1});
else
    bO  = f.bO(p_prop);
    muO = f.muO(p_prop);
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
vO     = -s.faceUpstr(upco, mobO).*s.T.*dpO;
bOvO   =  s.faceUpstr(upco, bO).*vO;
if disgas,
    rsbOvO = s.faceUpstr(upco, rs).*bOvO;
end

% GAS PROPS (calculated at oil pressure)
if vapoil
    bG  = f.bG(p_prop, rv, ~st{2});
    muG = f.muG(p_prop, rv, ~st{2});
else
    bG  = f.bG(p_prop);
    muG = f.muG(p_prop);
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
vG     = -s.faceUpstr(upcg, mobG).*s.T.*dpG;
bGvG   =  s.faceUpstr(upcg, bG).*vG;
if vapoil, 
    rvbGvG = s.faceUpstr(upcg, rv).*bGvG;
end

% These are needed in transport solver, so we output them regardless of
% any flags set in the model.
state = model.storeFluxes(state, vW, vO, vG);
state = model.storeUpstreamIndices(state, upcw, upco, upcg);
if model.extraStateOutput
    state = model.storebfactors(state, bW, bO, bG);
    state = model.storeMobilities(state, mobW, mobO, mobG);
end
% EQUATIONS -----------------------------------------------------------

bW0 = f.bW(p0);
if disgas, bO0 = f.bO(p0, rs0, ~st0{1}); else bO0 = f.bO(p0); end
if vapoil, bG0 = f.bG(p0, rv0, ~st0{2}); else bG0 = f.bG(p0); end

% oil eq:
if vapoil
    oil = (s.pv/dt).*( pvMult.* (bO.* sO  + rv.* bG.* sG) - ...
                          pvMult0.*(bO0.*sO0 + rv0.*bG0.*sG0) ) + ...
             s.Div(bOvO + rvbGvG);
else
    oil = (s.pv/dt).*( pvMult.*bO.*sO - pvMult0.*bO0.*sO0 ) + s.Div(bOvO);

end


% water eq:
wat = (s.pv/dt).*( pvMult.*bW.*sW - pvMult0.*bW0.*sW0 ) + s.Div(bWvW);

% gas eq:
if disgas
    gas = (s.pv/dt).*( pvMult.* (bG.* sG  + rs.* bO.* sO) - ...
                          pvMult0.*(bG0.*sG0 + rs0.*bO0.*sO0 ) ) + ...
             s.Div(bGvG + rsbOvO);
else
    gas = (s.pv/dt).*( pvMult.*bG.*sG - pvMult0.*bG0.*sG0 ) + s.Div(bGvG);
end



[eqs, names, types] = deal(cell(1, 5));
% well equations
if ~isempty(W)
    wm = WellModel();
    wc    = vertcat(W.cells);
    nperf = numel(wc);
    pw    = p(wc);
    rhows = [f.rhoWS, f.rhoOS, f.rhoGS];
    
    if ~disgas
        rsw = ones(nperf,1)*rs; rsSatw = ones(nperf,1)*rsSat; %constants
    else
        rsw = rs(wc); rsSatw = rsSat(wc);
    end
    if ~vapoil
        rvw = ones(nperf,1)*rv; rvSatw = ones(nperf,1)*rvSat; %constants
    else
        rvw = rv(wc); rvSatw = rvSat(wc);
    end
    rw    = {rsw, rvw};
    rSatw = {rsSatw, rvSatw};
    
    
    if 0
        sG_w = sG(wc) + dt*qGs./(bG(wc).*s.pv(wc));
        sO_w = sO(wc) + dt*qOs./(bO(wc).*s.pv(wc));
        sW_w = sW(wc) + dt*qWs./(bW(wc).*s.pv(wc));
        
        sG_w = min(max(sG_w, 0), 1);
        sO_w = min(max(sO_w, 0), 1);
        sW_w = 1 - sG_w - sO_w;
%         sW_w = min(max(sW_w, 0), 1);
        
        pw = p_prop(wc);
        
        
        sub = double(sG_w) > 0;
        rw{1}(sub) = f.rsSat(pw(sub));        
        
        if 1
            muG_w = muG(wc);
            muW_w = muW(wc);
            muO_w = muO(wc);
        else
            muG_w = f.muG(pw);
            if model.disgas
    %             muO_w = f.muO(pw, rw{1}, sG_w > 0);
                muO_w = f.muO(pw, rs(wc), ~st{1}(wc));
            else
                muO_w = f.muO(pw);
            end
            muW_w = f.muW(pw);
        end
        [krW_w, krO_w, krG_w] = model.evaluteRelPerm({sW_w, sO_w, sG_w});
        
        mobWw = krW_w./muW_w;
        mobOw = krO_w./muO_w;
        mobGw = krG_w./muG_w;
        
        mw    = {mobWw, mobOw, mobGw};
        bw    = {bW(wc), f.bO(pw, rw{1}, sG_w > 0), bG(wc)};
%         bw    = {bW(wc), bO(wc), bG(wc)};
    else
        mw    = {mobW(wc), mobO(wc), mobG(wc)};
        bw    = {bW(wc), bO(wc), bG(wc)};
    end
    
    s = {sW, 1 - sW - sG, sG};

    [cqs, weqs, ctrleqs, wc, state.wellSol, cqr]  = wm.computeWellFlux(model, W, wellSol, ...
                                         bhp, {qWs, qOs, qGs}, pw, rhows, bw, mw, s, rw,...
                                         'maxComponents', rSatw, ...
                                         'nonlinearIteration', opt.iteration);
    eqs(2:4) = weqs;
    eqs{5} = ctrleqs;

    qW = double(cqr{1});
    qO = double(cqr{2});
    qG = double(cqr{3});

    oil(wc) = oil(wc) - cqs{2}; % Add src to oil eq
    wat(wc) = wat(wc) - cqs{1}; % Add src to water eq
    gas(wc) = gas(wc) - cqs{3}; % Add src to gas eq

    names(2:5) = {'oilWells', 'waterWells', 'gasWells', 'closureWells'};
    types(2:5) = {'perf', 'perf', 'perf', 'well'};
end
% Create actual pressure equation
if disgas && vapoil
    cfac = 1./(1 - rs.*rv);
else
    cfac = 1;
end

if 1
    a_w = 1./bW;
    a_o = cfac.*(1./bO - rs./bG);
    a_g = cfac.*(1./bG - rv./bO);

    eqs{1} = oil.*a_o + wat.*a_w + gas.*a_g;
else
    a_w = 1./bW;
    a_o = 1./bO - rs./bG;
    a_g = 1./bG;

    wat = wat.*a_w;
    oil = oil.*a_o;
    gas = gas.*a_g;

    eqs{1} = wat + oil + gas;
end
names{1} = 'pressure';
types{1} = 'cell';

% Store fluxes for the transport solver
perf2well = getPerforationToWellMapping(W);
fluxt = qW + qO + qG;
for i = 1:numel(W)
    wp = perf2well == i;
    state.wellSol(i).flux = fluxt(wp);
end


% Hacking away...
state.s0 = state0.s;
state.bfactor0 = [double(bW0), double(bO0), double(bG0)];
% problem.state.state0 = state0;
% if isfield(problem.state.state0, 'state0')
%     problem.state.state0 = rmfield(problem.state.state0, 'state0');
% end

problem = LinearizedProblem(eqs, types, names, primaryVars, state, dt);

end


function [sG, rs] = redistributeRS(f, p, rs, sG, sO)
    rsSat = f.rsSat(p);
    isSat = rs >= rsSat;

    bG = f.bG(p);
    bO = f.bO(p, rs, isSat);

    assert(all(bO>0))

    % Find total Rs if everything was dissolved, i.e. sort of the mass
    % of gas for fixed compressibility
    dRs = sG.*bG./(max(double(sO), 0.001).*bO);
    rs = rs + dRs;
    rs(~isfinite(double(rs))) = 0;
    rs(double(rs)<0) = 0;

    sG = 0*sG;

    % Work out the overflow and put it into the gas phase
    above = rs>rsSat;
    overflow = rs(above) - rsSat(above);

    rs(above) = rsSat(above);

    sG(above) = overflow.*sO(above).*bO(above)./bG(above);
end

function s = estimateSaturationEndpoint()
    
end