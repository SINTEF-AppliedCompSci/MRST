function [problem, state] = pressureEquationBlackOil(state0, state, model, dt, drivingForces, varargin)

opt = struct('Verbose', mrstVerbose, ...
             'reverseMode', false,...
             'scaling', [],...
             'resOnly', false,...
             'history', [],...
             'iteration', -1, ...
             'stepOptions', []);  % Compatibility only

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
if ~opt.resOnly,
    if ~opt.reverseMode,
        % define primary varible x and initialize
        x = st{1}.*rs + st{2}.*rv + st{3}.*sG;

        [p, qWs, qOs, qGs, bhp] = ...
            initVariablesADI(p, qWs, qOs, qGs, bhp);
        % define sG, rs and rv in terms of x
        sG = st{2}.*(1-sW) + st{3}.*x;
        if disgas
            rsSat = f.rsSat(p);
            rs = (~st{1}).*rsSat + st{1}.*x;
        else % otherwise rs = rsSat = const
            rsSat = rs;
        end
        if vapoil
            rvSat = f.rvSat(p);
            rv = (~st{2}).*rvSat + st{2}.*x;
        else % otherwise rv = rvSat = const
            rvSat = rv;
        end
    else
        assert(0, 'Backwards solver not supported for splitting');
    end
else % resOnly-case compute rsSat and rvSat for use in well eqs
    if disgas, rsSat = f.rsSat(p); else rsSat = rs; end
    if vapoil, rvSat = f.rvSat(p); else rvSat = rv; end
end

primaryVars = {'pressure', 'qWs', 'qOs', 'qGs', 'bhp'};

    %----------------------------------------------------------------------
    %check for p-dependent tran mult:
    trMult = 1;
    if isfield(f, 'tranMultR'), trMult = f.tranMultR(p); end

    %check for p-dependent porv mult:
    pvMult = 1; pvMult0 = 1;
    if isfield(f, 'pvMultR')
        pvMult =  f.pvMultR(p);
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
    [krW, krO, krG] = f.relPerm(sW, sG);
    g  = norm(gravity);
    dz = s.grad(G.cells.centroids(:,3));

    % WATER PROPS (calculated at oil pressure)
    bW     = f.bW(p);
    rhoW   = bW.*f.rhoWS;
    % rhoW on face, avarge of neighboring cells (E100, not E300)
    rhoWf  = s.faceAvg(rhoW);
    mobW   = trMult.*krW./f.muW(p);
    dpW    = s.grad(p-pcOW) - g*(rhoWf.*dz);
    
    % water upstream-index
    upcw  = (double(dpW)>=0);
    vW   = s.faceUpstr(upcw, mobW).*s.T.*dpW;
    bWvW = s.faceUpstr(upcw, bW).*vW;
    if any(bW < 0)
        warning('Negative water compressibility present!')
    end
    
    % OIL PROPS
    if disgas
        bO  = f.bO(p, rs, ~st{1});
        muO = f.muO(p, rs, ~st{1});
    else
        bO  = f.bO(p);
        muO = f.muO(p);
    end
    if any(bO < 0)
        warning('Negative oil compressibility present!')
    end
    rhoO   = bO.*(rs*f.rhoGS + f.rhoOS);
    rhoOf  = s.faceAvg(rhoO);
    mobO   = trMult.*krO./muO;
    dpO    = s.grad(p) - g*(rhoOf.*dz);
    % oil upstream-index
    upco = (double(dpO)>=0);
    vO     = s.faceUpstr(upco, bO).*s.T.*dpO;
    bOvO   = s.faceUpstr(upco, mobO).*vO;
    if disgas,
        rsbOvO = s.faceUpstr(upco, rs).*bOvO;
    end

    % GAS PROPS (calculated at oil pressure)
    if vapoil
        bG  = f.bG(p, rv, ~st{2});
        muG = f.muG(p, rv, ~st{2});
    else
        bG  = f.bG(p);
        muG = f.muG(p);
    end
    if any(bG < 0)
        warning('Negative gas compressibility present!')
    end
    rhoG   = bG.*(rv*f.rhoOS + f.rhoGS);
    rhoGf  = s.faceAvg(rhoG);
    mobG   = trMult.*krG./muG;
    dpG    = s.grad(p+pcOG) - g*(rhoGf.*dz);
    % gas upstream-index
    upcg    = (double(dpG)>=0);
    vG     = s.faceUpstr(upcg, mobG).*s.T.*dpG;
    bGvG   = s.faceUpstr(upcg, bG).*vG;
    if vapoil, 
        rvbGvG = s.faceUpstr(upcg, rv).*bGvG;
    end

    % EQUATIONS -----------------------------------------------------------
    sO  = 1- sW  - sG;
    sO0 = 1- sW0 - sG0;

    bW0 = f.bW(p0);
    if disgas, bO0 = f.bO(p0, rs0, ~st0{1}); else bO0 = f.bO(p0); end
    if vapoil, bG0 = f.bG(p0, rv0, ~st0{2}); else bG0 = f.bG(p0); end

    % oil eq:
    if vapoil
        oil = (s.pv/dt).*( pvMult.* (bO.* sO  + rv.* bG.* sG) - ...
                              pvMult0.*(bO0.*sO0 + rv0.*bG0.*sG0) ) - ...
                 s.div(bOvO + rvbGvG);
    else
        oil = (s.pv/dt).*( pvMult.*bO.*sO - pvMult0.*bO0.*sO0 ) - s.div(bOvO);
        
    end
    
    
    % water eq:
    wat = (s.pv/dt).*( pvMult.*bW.*sW - pvMult0.*bW0.*sW0 ) - s.div(bWvW);
    
    % gas eq:
    if disgas
        gas = (s.pv/dt).*( pvMult.* (bG.* sG  + rs.* bO.* sO) - ...
                              pvMult0.*(bG0.*sG0 + rs0.*bO0.*sO0 ) ) - ...
                 s.div(bGvG + rsbOvO);
    else
        gas = (s.pv/dt).*( pvMult.*bG.*sG - pvMult0.*bG0.*sG0 ) - s.div(bGvG);
    end
    
    
    
    [eqs, names, types] = deal(cell(1, 5));
    % well equations
    if ~isempty(W)
        wm = WellModel();
        wc    = vertcat(W.cells);
        nperf = numel(wc);
        pw    = p(wc);
        rhows = [f.rhoWS, f.rhoOS, f.rhoGS];
        bw    = {bW(wc), bO(wc), bG(wc)};
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
        mw    = {mobW(wc), mobO(wc), mobG(wc)};
        s = {sW, 1 - sW - sG, sG};

        [cqs, weqs, ctrleqs, wc, state.wellSol]  = wm.computeWellFlux(model, W, wellSol, ...
                                             bhp, {qWs, qOs, qGs}, pw, rhows, bw, mw, s, rw,...
                                             'maxComponents', rSatw, ...
                                             'nonlinearIteration', opt.iteration);
        eqs(2:4) = weqs;
        eqs{5} = ctrleqs;

        
        qG = (cqs{3} - rsw.*cqs{2})./(1 - rsw.*rvw);
        qO = (cqs{2} - rvw.*cqs{3})./(1 - rsw.*rvw);
        
        qW = cqs{1};
%         qO = cqs{2};
%         qG = cqs{3};

        oil(wc) = oil(wc) - cqs{2}; % Add src to oil eq
        wat(wc) = wat(wc) - qW; % Add src to water eq
        gas(wc) = gas(wc) - cqs{3}; % Add src to gas eq

        names(2:5) = {'oilWells', 'waterWells', 'gasWells', 'closureWells'};
        types(2:5) = {'perf', 'perf', 'perf', 'well'};
    end
    % Create actual pressure equation
    cfac = 1./(1 - rs.*rv);
    a_o = cfac.*(1./bO - rs./bG);
    a_g = cfac.*(1./bG - rv./bO);
    a_w = 1./bW;
    
    eqs{1} = oil.*a_o + wat.*a_w + gas.*a_g;
    names{1} = 'pressure';
    types{1} = 'cell';
    
    % Store fluxes for the transport solver
    perf2well = getPerforationToWellMapping(W);
    for i = 1:numel(W)
        wp = perf2well == i;
        state.wellSol(i).flux = double(qW(wp)) + double(qO(wp)) + double(qG(wp));
    end
    

    
    state.upstream = [upcw, upco, upcg];
    
    intx = ~any(G.faces.neighbors == 0, 2);
    if ~isfield(state, 'flux')
        state.flux = zeros(G.faces.num, 1);
    end
    state.flux(intx) = double(vW + vO + vG);
    
    problem = LinearizedProblem(eqs, types, names, primaryVars, state, dt);
    problem.iterationNo = opt.iteration;
end
