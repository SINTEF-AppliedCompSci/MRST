function [problem, state] = transportEquationOilWaterDG(state0, state, model, dt, drivingForces, varargin)

    opt = struct('Verbose', mrstVerbose, ...
                 'reverseMode', false,...
                 'scaling', [],...
                 'resOnly', false,...
                 'history', [],...
                 'solveForWater', false, ...
                 'solveForOil', true, ...
                 'iteration', -1, ...
                 'stepOptions', []);  % Compatibility only
    
    opt   = merge_options(opt, varargin{:});
    W     = drivingForces.W;
    op    = model.operators;
    fluid = model.fluid;
    rock  = model.rock;
    G     = model.G;
    disc  = model.disc;
        
    assert(~(opt.solveForWater && opt.solveForOil));
    
    
    if opt.iteration == 1
        
        state.degree = repmat(disc.degree, G.cells.num, 1);
        state        = disc.mapDofs(state, state0);
        
%         state0_tmp = state0;
%         state0.degree = repmat(disc.degree, G.cells.num, 1);
%         state0        = disc.mapDofs(state0, state0_tmp);
        
        
    end
    state = disc.updateDisc(state);
    

    [p , sWdof , wellSol] = model.getProps(state , 'pressure', 'water', 'wellsol');
    [p0, sWdof0         ] = model.getProps(state0, 'pressure', 'water'           );

    % If timestep has been split relative to pressure, linearly interpolate in
    % pressure.
    pFlow = p;
    if isfield(state, 'timestep')
        dt_frac = dt/state.timestep;
        p = p.*dt_frac + p0.*(1-dt_frac);
    end
    
    %Initialization of independent variables ------------------------------

    if ~opt.resOnly,
        % ADI variables needed since we are not only computing residuals.
        if ~opt.reverseMode,
            sWdof = model.AutoDiffBackend.initVariablesAD(sWdof);
        else
            assert(0, 'Backwards solver not supported for splitting');
        end
    end
    
    % ---------------------------------------------------------------------

    primaryVars = {'sWdof'};
    
%     nDof    = disc.getnDof();

%     nDof = disc.nDof;
    nDof = state.nDof;
    
    % Express sW and sW0 in basis
    sW  = @(x,c) disc.evaluateSaturation(x, c, sWdof , state );
    sW0 = @(x,c) disc.evaluateSaturation(x, c, sWdof0, state0);
    sO  = @(x,c) 1-sW(x,c);
    
    [pvMult, transMult, mobMult, pvMult0] = getMultipliers(model.fluid, p, p0);
    T = op.T.*transMult;
    gdz = 0;
    [vW, bW, mobW, rhoW, pW, upcW, dpW, muW] = getPropsWater_DG(model, p, sW, T, gdz);
    bW0 = fluid.bW(p0);
    
    [vO, bO, mobO, rhoO, pO, upcO, dpO, muO] = getPropsOil_DG(model, p, sO, T, gdz);
    
    % Accumulation term----------------------------------------------------
    
    if numel(pvMult) == 1
        pvMult = repmat(pvMult, G.cells.num,1);
    end
    if numel(pvMult0) == 1
        pvMult0 = repmat(pvMult0, G.cells.num,1);
    end
    
    
%     acc = @(sW, sW0, c, psi) (pvMult (c).*rock.poro(c).*bW (c).*sW - ...
%                               pvMult0(c).*rock.poro(c).*bW0(c).*sW0).*psi/dt;
    
    acc = @(x,c,psi) (pvMult (c).*rock.poro(c).*bW (c).*sW (x,c) - ...
                            pvMult0(c).*rock.poro(c).*bW0(c).*sW0(x,c)).*psi/dt;
    
    % Flux term------------------------------------------------------------
    
    vT  = sum(state.flux,2);
    vTc = faceFlux2cellVelocity(G, vT);
    
    gp = op.Grad(p);
    
    [Gw, Go] = deal(zeros(G.faces.num, 1));
    Gw(op.internalConn) = op.T.*(gp - dpW);
    Go(op.internalConn) = op.T.*(gp - dpO);
    
    Gwc = faceFlux2cellVelocity(G, Gw);
    Goc = faceFlux2cellVelocity(G, Go);

    fW = @(x,c) mobW(x,c)./(mobW(x,c) + mobO(x,c));

    flux1 = @(x,c,grad_psi) bW(c).*fW(x, c).*sum(vTc(c,:).*grad_psi,2) ...
                 + bO(c).*fW(x, c).*sum((Gwc(c,:) - Goc(c,:)).*grad_psi,2);
    
%     flux1 = @(x,c,grad_psi) bW(c).*fW(x, c).*sum(vTc(c,:).*grad_psi,2) ...
%                               + bO(c).*fW(x, c).*sum((Gwc(c,:) - Goc(c,:)).*grad_psi,2);
                          
    cellIntegral = disc.cellInt(@(x,c,psi,grad_psi) acc(x,c,psi) - flux1(x,c,grad_psi), (1:G.cells.num)', state);
    
%     flux1 = -disc.cellIntDiv(flux1, (1:G.cells.num)', state);
    
    
    
    flux2 = @(xc, xv, xg, c, cv, cg, f, psi) ...
        (bW(cg).*fW(xv, cv).*vT(f)./G.faces.areas(f) ...
       + bO(cg).*fW(xg, cg).*mobO(xg,cg).*(Gw(f) - Go(f))).*psi;

%     flux2 = disc.faceIntDiv(integrand, (1:G.cells.num)', upcW, state);
    faceIntegral = disc.faceIntDiv(flux2, (1:G.cells.num)', upcW, state);
  
    % Water equation-------------------------------------------------------

    water = cellIntegral + faceIntegral;

    % Well contributions---------------------------------------------------
    
    if ~isempty(W)
        
        perf2well = getPerforationToWellMapping(W);
        wc = vertcat(W.cells);

        wflux = zeros(G.cells.num,1);
        wflux(wc) = sum(vertcat(wellSol.flux), 2);
        isInj = wflux > 0;
        compWell = vertcat(W.compi);
        compPerf = zeros(G.cells.num, 2);
        compPerf(wc,:) = compWell(perf2well,:);
        
        integrand = @(x, c, psi, grad_psi) ...
            bW(c).*wflux(c).*(fW(x, c)     .*(~isInj(c)) ...
                            + compPerf(c,1).*( isInj(c))).*psi;
        
        prod = disc.cellInt(integrand, wc, state);
        
%         vol = reshape(repmat(G.cells.volumes(wc)', nDofMax, 1), [], 1);
        
        vol = rldecode(G.cells.volumes(wc), nDof(wc), 1);
        
        ix = disc.getDofIx(state, [], wc);
        water(ix) = water(ix) - prod(ix)./vol;

        % Store well fluxes
        [x, w, nq, ii, jj, cellNo] = makeCellIntegrator(G, wc, disc.degree+1, 'volume');
        x = disc.transformCoords(x, cellNo);
        WC = sparse(ii, jj, w);
        
        wflux_W = WC*(bW(cellNo).*wflux(cellNo) ...
                     .*(fW(x, cellNo) .*(~isInj(cellNo)) ...
                     + compPerf(cellNo,1).*( isInj(cellNo))));
        wflux_W = wflux_W./G.cells.volumes(wc);
          
        wflux_O = WC*(bO(cellNo).*wflux(cellNo) ...
                     .*((1-fW(x, cellNo)).*(~isInj(cellNo)) ...
                     + compPerf(cellNo,2) .*( isInj(cellNo))));
        wflux_O = wflux_O./G.cells.volumes(wc);

        wflux_O = double(wflux_O);
        wflux_W = double(wflux_W);

        for wNo = 1:numel(W)
            perfind = perf2well == wNo;
            state.wellSol(wNo).qOs = sum(wflux_O(perfind));
            state.wellSol(wNo).qWs = sum(wflux_W(perfind));
        end

    end
    
    state.cfl = dt.*sum(abs(vTc)./G.cells.dx,2);
    
    % Make Linearized problem----------------------------------------------

    eqs   = {water  };
    names = {'water'};
    types = {'cell' };

    if ~model.useCNVConvergence
        pv = rldecode(op.pv, nDof, 1);
        eqs{1} = eqs{1}.*(dt./pv);
    end

    problem = LinearizedProblem(eqs, types, names, primaryVars, state, dt);

end