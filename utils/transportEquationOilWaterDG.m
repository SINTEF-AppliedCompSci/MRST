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

    [p , sWdof , wellSol] = model.getProps(state , 'pressure', 'water', 'wellsol');
    [p0, sWdof0         ] = model.getProps(state0, 'pressure', 'water'           );

    % If timestep has been split relative to pressure, linearly interpolate in
    % pressure.
    pFlow = p;
    if isfield(state, 'timestep')
        dt_frac = dt/state.timestep;
        p = p.*dt_frac + p0.*(1-dt_frac);
    end
    
    %Initialization of independent variables ----------------------------------

    if ~opt.resOnly,
        % ADI variables needed since we are not only computing residuals.
        if ~opt.reverseMode,
            sWdof = model.AutoDiffBackend.initVariablesAD(sWdof);
        else
            assert(0, 'Backwards solver not supported for splitting');
        end
    end
    
    % -------------------------------------------------------------------------

    primaryVars = {'sWdof'};
    
%     [psi, grad_psi, k, nDof] = dgBasis(model.degree, model.G.griddim, 'legendre');

    psi      = disc.basis.psi;
    grad_psi = disc.basis.grad_psi;
    nDof     = disc.basis.nDof;
    
    % Express sW and sW0 in basis
    sW  = @(x,c) getSatFromDof(x, c, sWdof , disc);
    sW0 = @(x,c) getSatFromDof(x, c, sWdof0, disc);
    sO  = @(x,c) 1-sW(x,c);
    
    [pvMult, transMult, mobMult, pvMult0] = getMultipliers(model.fluid, p, p0);
    T = op.T.*transMult;
    gdz = 0;
    [vW, bW, mobW, rhoW, pW, upcW, dpW, muW] = getPropsWater_DG(model, p, sW, T, gdz);
    bW0 = fluid.bW(p0);
    
    [vO, bO, mobO, rhoO, pO, upcO, dpO, muO] = getPropsOil_DG(model, p, sO, T, gdz);
    
    [xc, cellNo_c,         WC] = cellBasisIntegrator(disc);
    [xf, cellNo_f, faceNo, WF] = faceBasisIntegrator(disc);
    
    % Accumulation term----------------------------------------------------
    
    if numel(pvMult) == 1
        pvMult = repmat(pvMult, G.cells.num,1);
    end
    if numel(pvMult0) == 1
        pvMult0 = repmat(pvMult0, G.cells.num,1);
    end
    
    acc  = sWdof;
    for dofNo = 1:nDof
        
        ix = (1:nDof:G.cells.num*nDof) + dofNo - 1;
        now  = WC*(pvMult(cellNo_c) .*bW(cellNo_c) .*rock.poro(cellNo_c).*sW(xc,cellNo_c) .*psi{dofNo}(xc));
        then = WC*(pvMult0(cellNo_c).*bW0(cellNo_c).*rock.poro(cellNo_c).*sW0(xc,cellNo_c).*psi{dofNo}(xc));
        acc(ix) = (now - then)/dt;
        
    end
    
    % Flux term------------------------------------------------------------
    
    vT = sum(state.flux,2);
    vTc = faceFlux2cellVelocity(G, vT);
    
    gp = op.Grad(p);
    
    [Gw, Go] = deal(zeros(G.faces.num, 1));
    Gw(op.internalConn) = op.T.*(gp - dpW);
    Go(op.internalConn) = op.T.*(gp - dpO);
    
    Gwc = faceFlux2cellVelocity(G, Gw);
    Goc = faceFlux2cellVelocity(G, Go);

    fW = @(x,c) mobW(x,c)./(mobW(x,c) + mobO(x,c));
    
    flux1 = sWdof;
    for dofNo = 1:nDof
        
        ix = (1:nDof:G.cells.num*nDof) + dofNo - 1;
        flux1(ix) = -WC*(bW(cellNo_c).*fW(xc, cellNo_c).*sum(vTc(cellNo_c,:).*grad_psi{dofNo}(xc),2)  ...
                       + bO(cellNo_c).*fW(xc, cellNo_c).*sum((Gwc(cellNo_c,:) - Goc(cellNo_c,:)).*grad_psi{dofNo}(xc),2));
                   
    end
    
    upCells_v = G.faces.neighbors(:,2);
    
    intf = find(op.internalConn);
    upCells_v(intf(upcW)) = op.N(upcW,1);

    upCells_v = upCells_v(faceNo);    
    upCells_G = upCells_v;
    
    flux2 = sWdof;
    for dofNo = 1:nDof
        
        ix = (1:nDof:G.cells.num*nDof) + dofNo - 1;
        flux2(ix) = WF*(bW(upCells_G).*fW(xf, upCells_v).*vT(faceNo).*psi{dofNo}(xf) ...
                      + bO(upCells_G).*fW(xf, upCells_G).*mobO(xf,upCells_G).*(Gw(faceNo) - Go(faceNo)).*psi{dofNo}(xf));
                  
    end
    
    flux  = flux1 + flux2;
    water = acc + flux;
    
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

        [ii, jj] = find(WC);
        keep = any(ii == wc',2);
        jj = jj(keep);
        
        S1 = sparse((1:numel(wc))', wc, 1, numel(wc), G.cells.num);
        S2 = sparse(jj, (1:numel(jj))' , 1, size(WC,2)     , numel(jj));
        WWC = S1*WC*S2;
        
        keep = any(cellNo_c == wc',2);
        xwc = xc(keep,:);
        cellNo_wc = cellNo_c(keep);
        
        prod = sWdof(wc);
        for dofNo = 1:nDof
            ix = (1:nDof:numel(wc)*nDof) + dofNo - 1;
            prod(ix) = (WWC*(bW(cellNo_wc).*wflux(cellNo_wc)...
                          .*(fW(xwc, cellNo_wc) .*(~isInj(cellNo_wc)) ...
                          +  compPerf(cellNo_wc,1).*( isInj(cellNo_wc))).*psi{dofNo}(xwc)))./G.cells.volumes(wc);
        end
        
        ind = mcolon((wc-1)*nDof + 1, wc*nDof);
        water(ind) = water(ind) - prod;

        % Store well fluxes
        wflux_W = bW(wc).*wflux(wc) ...
              .*(fW([0,0], wc) .*(~isInj(wc)) ...
               + compPerf(wc,1).*( isInj(wc)));
          
        wflux_O = bO(wc).*wflux(wc) ...
              .*((1-fW([0,0], wc)).*(~isInj(wc)) ...
               + compPerf(wc,2) .*( isInj(wc)));
          
          
        
        wflux_O = double(wflux_O);
        wflux_W = double(wflux_W);
% 
        for wNo = 1:numel(W)
            perfind = perf2well == wNo;
            state.wellSol(wNo).qOs = sum(wflux_O(perfind));
            state.wellSol(wNo).qWs = sum(wflux_W(perfind));
        end

    end

eqs = {water};
names = {'water'};
types = {'cell'};

pv = reshape(repmat(op.pv', nDof, 1), [], 1);
if ~model.useCNVConvergence
    eqs{1} = eqs{1}.*(dt./pv);
end

problem = LinearizedProblem(eqs, types, names, primaryVars, state, dt);

if 1
    faces = G.cells.faces(:,1);
    cells = rldecode((1:G.cells.num)', diff(G.cells.facePos), 1);
    xx = (G.faces.centroids(faces,:) - G.cells.centroids(cells,:))./(2*sqrt(G.griddim)*G.cells.diameters(cells,:));
    ss = sW(xx,cells);
end

end