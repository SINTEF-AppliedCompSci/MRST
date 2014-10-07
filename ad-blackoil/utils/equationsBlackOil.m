function [problem, state] = equationsBlackOil(state0, state, model, dt, drivingForces, varargin)
% Generate linearized problem for a volatile 3Ph system (wet-gas, live-oil).
    opt = struct('Verbose',     mrstVerbose,...
                 'reverseMode', false,...
                 'resOnly',     false,...
                 'iteration',   -1);

    opt = merge_options(opt, varargin{:});
    
    W = drivingForces.Wells;
    assert(isempty(drivingForces.bc) && isempty(drivingForces.src))
    
    s = model.operators;
    G = model.G;
    f = model.fluid;

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
    
    %Initialization of primary variables ----------------------------------
    st  = getCellStatusVO(state,  1-sW-sG,   sW,  sG,  disgas, vapoil);
    st0 = getCellStatusVO(state0, 1-sW0-sG0, sW0, sG0, disgas, vapoil);
    if ~opt.resOnly,
        if ~opt.reverseMode,
            % define primary varible x and initialize
            x = st{1}.*rs + st{2}.*rv + st{3}.*sG;

            [p, sW, x, qWs, qOs, qGs, bhp] = ...
                initVariablesADI(p, sW, x, qWs, qOs, qGs, bhp);
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
            x0 = st0{1}.*rs0 + st0{2}.*rv0 + st0{3}.*sG0;

            [p0, sW0, x0, zw, zw, zw, zw] = ...
                initVariablesADI(p0, sW0, x0, ...
                zeros(size(qWs)) , zeros(size(qOs)) , ...
                zeros(size(qGs)) , zeros(size(bhp)));                %#ok
            sG0 = st0{2}.*(1-sW0) + st0{3}.*x0;
            if disgas
                rsSat0 = f.rsSat(p0);
                rs0 = (~st0{1}).*rsSat0  + st0{1}.*x0;
            end
            if vapoil
                rvSat0 = f.rvSat(p0);
                rv0 = (~st0{2}).*rvSat0  + st0{2}.*x0;
            end
        end
    else % resOnly-case compute rsSat and rvSat for use in well eqs
        if disgas, rsSat = f.rsSat(p); else rsSat = rs; end
        if vapoil, rvSat = f.rvSat(p); else rvSat = rv; end
    end
    if disgas || vapoil
        gvar = 'x';
    else
        gvar = 'sG';
    end
    primaryVars = {'pressure', 'sW', gvar, 'qWs', 'qOs', 'qGs', 'bhp'};
    
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
    vW = s.faceUpstr(upcw, mobW).*s.T.*dpW;
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
    vO   = s.faceUpstr(upco, mobO).*s.T.*dpO;
    bOvO   = s.faceUpstr(upco, bO).*vO;
    if disgas
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
    vG = s.faceUpstr(upcg, mobG).*s.T.*dpG;
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
    sO  = 1- sW  - sG;
    sO0 = 1- sW0 - sG0;

    bW0 = f.bW(p0);
    if disgas, bO0 = f.bO(p0, rs0, ~st0{1}); else bO0 = f.bO(p0); end
    if vapoil, bG0 = f.bG(p0, rv0, ~st0{2}); else bG0 = f.bG(p0); end

    % oil eq:
    names{1} = 'oil';
    if vapoil
        eqs{1} = (s.pv/dt).*( pvMult.* (bO.* sO  + rv.* bG.* sG) - ...
                              pvMult0.*(bO0.*sO0 + rv0.*bG0.*sG0) ) - ...
                 s.div(bOvO + rvbGvG);
    else
        eqs{1} = (s.pv/dt).*( pvMult.*bO.*sO - pvMult0.*bO0.*sO0 ) - s.div(bOvO);
        
    end
    
    
    % water eq:
    names{2} = 'water';
    eqs{2} = (s.pv/dt).*( pvMult.*bW.*sW - pvMult0.*bW0.*sW0 ) - s.div(bWvW);
    
    % gas eq:
    names{3} = 'gas';
    if disgas
        eqs{3} = (s.pv/dt).*( pvMult.* (bG.* sG  + rs.* bO.* sO) - ...
                              pvMult0.*(bG0.*sG0 + rs0.*bO0.*sO0 ) ) - ...
                 s.div(bGvG + rsbOvO);
    else
        eqs{3} = (s.pv/dt).*( pvMult.*bG.*sG - pvMult0.*bG0.*sG0 ) - s.div(bGvG);
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
            eqs(4:6) = weqs;
            eqs{7} = ctrleqs;
            
            eqs{1}(wc) = eqs{1}(wc) - cqs{2}; % Add src to oil eq
            eqs{2}(wc) = eqs{2}(wc) - cqs{1}; % Add src to water eq
            eqs{3}(wc) = eqs{3}(wc) - cqs{3}; % Add src to gas eq
            
            names(4:7) = {'oilWells', 'waterWells', 'gasWells', 'closureWells'};
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
    problem.iterationNo = opt.iteration;
end

