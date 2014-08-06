function [problem, state] = equationsBlackOil(state0, state, model, dt, drivingForces, varargin)
% Generate linearized problem for a volatile 3Ph system (wet-gas, live-oil).
    opt = struct('Verbose',     mrstVerbose,...
                 'reverseMode', false,...
                 'scaling',     [],...
                 'resOnly',     false,...
                 'history',     [],  ...
                 'iteration',   -1,  ...
                 'stepOptions', []);

    opt = merge_options(opt, varargin{:});
    
    W = drivingForces.Wells;
    assert(isempty(drivingForces.bc) && isempty(drivingForces.src))
    
    s = model.operators;
    G = model.G;
    f = model.fluid;

    disgas = model.disgas;
    vapoil = model.vapoil;

    % Oil pressure
    p  = model.getProp(state, 'pressure');
    p0 = model.getProp(state0, 'pressure');

    % Water saturation
    sW  = model.getProp(state,  'water');
    sW0 = model.getProp(state0, 'water');
    
    % Gas saturation
    sG  = model.getProp(state,  'gas');
    sG0 = model.getProp(state0, 'gas');
    
    % Gas component in oil phase
    rs  = model.getProp(state,  'rs');
    rs0 = model.getProp(state0, 'rs');
    
    % Oil component in gas phase
    rv  = model.getProp(state,  'rv');
    rv0 = model.getProp(state0, 'rv');


    bhp = vertcat(state.wellSol.bhp);
    qWs    = vertcat(state.wellSol.qWs);
    qOs    = vertcat(state.wellSol.qOs);
    qGs    = vertcat(state.wellSol.qGs);
    
    %Initialization of primary variables ----------------------------------
    [st1 , st2  , st3 ] = getCellStatus(state , disgas, vapoil);
    [st1p, st2p , st3p] = getCellStatus(state0, disgas, vapoil);
    if ~opt.resOnly,
        if ~opt.reverseMode,
            % define primary varible x and initialize
            x = st1.*rs + st2.*rv + st3.*sG;

            [p, sW, x, qWs, qOs, qGs, bhp] = ...
                initVariablesADI(p, sW, x, qWs, qOs, qGs, bhp);
            % define sG, rs and rv in terms of x
            sG = st2.*(1-sW) + st3.*x;
            if disgas
                rsSat = f.rsSat(p);
                rs = (~st1).*rsSat + st1.*x;
            else % otherwise rs = rsSat = const
                rsSat = rs;
            end
            if vapoil
                rvSat = f.rvSat(p);
                rv = (~st2).*rvSat + st2.*x;
            else % otherwise rv = rvSat = const
                rvSat = rv;
            end
        else
            x0 = st1p.*rs0 + st2p.*rv0 + st3p.*sG0;

            [p0, sW0, x0, zw, zw, zw, zw] = ...
                initVariablesADI(p0, sW0, x0, ...
                zeros(size(qWs)) , zeros(size(qOs)) , ...
                zeros(size(qGs)) , zeros(size(bhp)));                %#ok
            sG0 = st2p.*(1-sW0) + st3p.*x0;
            if disgas
                rsSat0 = f.rsSat(p0);
                rs0 = (~st1p).*rsSat0  + st1p.*x0;
            end
            if vapoil
                rvSat0 = f.rvSat(p0);
                rv0 = (~st2p).*rvSat0  + st2p.*x0;
            end
        end
    else % resOnly-case compute rsSat and rvSat for use in well eqs
        if disgas, rsSat = f.rsSat(p); else rsSat = rs; end
        if vapoil, rvSat = f.rvSat(p); else rvSat = rv; end
    end
    primaryVars = {'pressure', 'sW', 'x', 'qWs', 'qOs', 'qGs', 'bhp'};
    
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
    upc  = (double(dpW)>=0);
    bWvW = s.faceUpstr(upc, bW.*mobW).*s.T.*dpW;
    if any(bW < 0)
        warning('Negative water compressibility present!')
    end
    
    % OIL PROPS
    if disgas
        bO  = f.bO(p, rs, ~st1);
        muO = f.muO(p, rs, ~st1);
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
    upc = (double(dpO)>=0);
    bOvO   = s.faceUpstr(upc, bO.*mobO).*s.T.*dpO;
    if disgas, rsbOvO = s.faceUpstr(upc, rs).*bOvO;end

    % GAS PROPS (calculated at oil pressure)
    if vapoil
        bG  = f.bG(p, rv, ~st2);
        muG = f.muG(p, rv, ~st2);
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
    upc    = (double(dpG)>=0);
    bGvG   = s.faceUpstr(upc, bG.*mobG).*s.T.*dpG;
    if vapoil, rvbGvG = s.faceUpstr(upc, rv).*bGvG; end

    % EQUATIONS -----------------------------------------------------------
    sO  = 1- sW  - sG;
    sO0 = 1- sW0 - sG0;

    bW0 = f.bW(p0);
    if disgas, bO0 = f.bO(p0, rs0, ~st1p); else bO0 = f.bO(p0); end
    if vapoil, bG0 = f.bG(p0, rv0, ~st2p); else bG0 = f.bG(p0); end

    % oil eq:
    names{1} = 'oil';
    if vapoil
        eqs{1} = (s.pv/dt).*( pvMult.* (bO.* sO  + rv.* bG.* sG) - ...
                              pvMult0.*(bO0.*sO0 + rv0.*bG0.*sG0) ) + ...
                 s.div(bOvO + rvbGvG);
    else
        eqs{1} = (s.pv/dt).*( pvMult.*bO.*sO - pvMult0.*bO0.*sO0 ) + s.div(bOvO);
        
    end
    
    
    % water eq:
    names{2} = 'water';
    eqs{2} = (s.pv/dt).*( pvMult.*bW.*sW - pvMult0.*bW0.*sW0 ) + s.div(bWvW);
    
    % gas eq:
    names{3} = 'gas';
    if disgas
        eqs{3} = (s.pv/dt).*( pvMult.* (bG.* sG  + rs.* bO.* sO) - ...
                              pvMult0.*(bG0.*sG0 + rs0.*bO0.*sO0 ) ) + ...
                 s.div(bGvG + rsbOvO);
    else
        eqs{3} = (s.pv/dt).*( pvMult.*bG.*sG - pvMult0.*bG0.*sG0 ) + s.div(bGvG);
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
            
            [cqs, weqs, ctrleqs, wc, state.wellSol]  = wm.computeWellFlux(model, W, state.wellSol, ...
                                                 bhp, {qWs, qOs, qGs}, pw, rhows, bw, mw, s,...
                                                 'pseudocomponents',    rw, ...
                                                 'maxPseudocomponents', rSatw, ...
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
    else
        eqs(4:7) = {bhp, bhp, bhp, bhp};  % empty  ADIs
        names(4:7) = {'empty', 'empty', 'empty', 'empty'};
        types(4:7) = {'none', 'none', 'none', 'none'};
    end
     
    
    problem = LinearizedProblem(eqs, types, names, primaryVars, state, dt);
    problem.iterationNo = opt.iteration;
end
%--------------------------------------------------------------------------
function [st1, st2, st3] = getCellStatus(state, disgas, vapoil)
% Status should be passed on from updateStateVO (to be sure definition is
% identical). rs and rv are assumed to be compatible, i.e. rx = rxSat for
% saturated cells and rx <= rxSat for undersaturated. Three values of
% status are:
% status 0: should not occur (almost water only -> state 3)
% status 1 oil, no gas  : x = rs, sg = 0    , rv = rvMax
% status 2 gas, no oil  : x = rv, sg = 1-sw , rs = rsMax
% status 3 oil and gas  : x = sg, rs = rsMax, rv = rvMax
if isfield(state, 'status')
    status = state.status;
else
    s = state.s;
    watOnly    = s(:,1) > 1- sqrt(eps);
    if ~vapoil
        oilPresent = true;
    else
        oilPresent = or(s(:,2) > 0, watOnly);
    end
    if ~disgas
        gasPresent = true;
    else
        gasPresent = or(s(:,3) > 0, watOnly);
    end
    status = oilPresent + 2*gasPresent;
end
if ~disgas
    st1 = false;
else
    st1 = status==1;
end
if ~vapoil
    st2 = false;
else
    st2 = status==2;
end
st3 = status == 3;
end

