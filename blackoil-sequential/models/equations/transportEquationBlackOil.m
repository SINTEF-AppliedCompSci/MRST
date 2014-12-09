function [problem, state] = transportEquationBlackOil(state0, state, model, dt, drivingForces, varargin)
% Generate linearized problem for a volatile 3Ph system (wet-gas, live-oil).
    opt = struct('Verbose',     mrstVerbose,...
                 'reverseMode', false,...
                 'scaling',     [],...
                 'resOnly',     false,...
                 'history',     [],  ...
                 'solveForWater', false, ...
                 'solveForOil', true, ...
                 'solveForGas', true, ...
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

    % Properties at current timestep
    [p, sW, sG, rs, rv, wellSol] = model.getProps(state, ...
                                    'pressure', 'water', 'gas', 'rs', 'rv', 'wellSol');
    % Properties at previous timestep
    [p0, sW0, sG0, rs0, rv0] = model.getProps(state0, ...
                                    'pressure', 'water', 'gas', 'rs', 'rv');

    wflux = vertcat(wellSol.flux);

    %Initialization of primary variables ----------------------------------
    st  = getCellStatusVO(state,  1-sW-sG,   sW,  sG,  disgas, vapoil);
    st0 = getCellStatusVO(state0, 1-sW0-sG0, sW0, sG0, disgas, vapoil);
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
    sO  = 1 - sW  - sG;
    sO0 = 1 - sW0 - sG0;

    [krW, krO, krG] = model.evaluteRelPerm({sW, sO, sG});
    % Gravity contribution
    gdz = model.getGravityGradient();
    
    % WATER PROPS (calculated at oil pressure)
    bW     = f.bW(p);
    rhoW   = bW.*f.rhoWS;
    % rhoW on face, avarge of neighboring cells (E100, not E300)
    rhoWf  = s.faceAvg(rhoW);
    mobW   = trMult.*krW./f.muW(p);
    

    
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


    % EQUATIONS -----------------------------------------------------------
    bW0 = f.bW(p0);
    if disgas, bO0 = f.bO(p0, rs0, ~st0{1}); else bO0 = f.bO(p0); end
    if vapoil, bG0 = f.bG(p0, rv0, ~st0{2}); else bG0 = f.bG(p0); end

    % Get total flux from state
    flux = sum(state.flux, 2);
    vT = flux(model.operators.internalConn);
    
    % Sat dependent pressure terms
    Go = rhoOf.*gdz;
    Gw = rhoWf.*gdz;
    Gg = rhoGf.*gdz;

    if numel(double(pcOW)) > 1
        Gw = Gw - s.Grad(pcOW);
    end
    
    if numel(double(pcOG)) > 1
        Gg = Gg - s.Grad(pcOG);
    end
    
    % Stored upstream indices
    if model.staticUpwind
        flag = state.upstreamFlag;
    else
        flag = multiphaseUpwindIndices({Gw, Go, Gg}, vT, s.T, ...
                {mobW, mobO, mobG}, s.faceUpstr);
    end

    upcw  = flag(:, 1);
    upco  = flag(:, 2);
    upcg  = flag(:, 3);
    
    % Upstream weighted face mobilities
    mobWf = s.faceUpstr(upcw, mobW);
    mobOf = s.faceUpstr(upco, mobO);
    mobGf = s.faceUpstr(upcg, mobG);
    % Tot mob
    totMob = mobOf + mobWf + mobGf;
    totMob = max(totMob, sqrt(eps));
    
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
    eqInd = 1;
    if opt.solveForWater
        % water eq:
        wat = (s.pv/dt).*( pvMult.*bW.*sW - pvMult0.*bW0.*sW0 ) + s.Div(bWvW);
        wat(wc) = wat(wc) - wflux_W;
        eqs{eqInd}   = wat;
        
        names{eqInd} = 'water';
        types{eqInd} = 'cell';
        eqInd = eqInd + 1;
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
        
        gas(wc) = gas(wc) - wflux_G;
        
        eqs{eqInd}   = gas;
        names{eqInd} = 'gas';
        types{eqInd} = 'cell';
        eqInd = eqInd + 1;
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
        oil(wc) = oil(wc) - wflux_O;
        
        eqs{eqInd}   = oil;
        names{eqInd} = 'oil';
        types{eqInd} = 'cell';
    end
    
    problem = LinearizedProblem(eqs, types, names, primaryVars, state, dt);
end

