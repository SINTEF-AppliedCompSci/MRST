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
    
    opt = merge_options(opt, varargin{:});
    W   = drivingForces.W;
    op   = model.operators;
    fluid = model.fluid;
    rock = model.rock;
    G = model.G;
        
    assert(~(opt.solveForWater && opt.solveForOil));

    [p, sWdof, wellSol] = model.getProps(state, 'pressure', 'water', 'wellsol');

    [p0, sWdof0] = model.getProps(state0, 'pressure', 'water');

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

    psi      = model.basis.psi;
    grad_psi = model.basis.grad_psi;
    nDof     = model.basis.nDof;
    
    % Express sW and sW0 in basis
    sW  = @(x,c) getSatFromDof(x, c, sWdof , model);
    sW0 = @(x,c) getSatFromDof(x, c, sWdof0, model);
    sO  = @(x,c) 1-sW(x,c);
    
    [pvMult, transMult, mobMult, pvMult0] = getMultipliers(model.fluid, p, p0);
    T = op.T.*transMult;
    gdz = 0;
    [vW, bW, mobW, rhoW, pW, upcW, dpW, muW] = getPropsWater_DG(model, p, sW, T, gdz);
    bW0 = fluid.bW(p0);
    
    [vO, bO, mobO, rhoO, pO, upcO, dpO, muO] = getPropsOil_DG(model, p, sO, T, gdz);
    
    [xc, wc, nqc, ii, jj, cellNoc] = makeCellIntegrator(G, (1:G.cells.num)', model.degree, 'tri');
    xc = (xc - G.cells.centroids(cellNoc,:))./G.cells.diameters(cellNoc);
    
    [xf, wf, cellsf, faces, nqf] = makeFaceIntegrator(G, (1:G.cells.num)', model.degree);
    xf = (xf - G.cells.centroids(cellsf,:))./G.cells.diameters(cellsf);
    
    wc = repmat(wc(:), nDof, 1);
    [ii, jj] = blockDiagIndex(ones(G.cells.num*nDof, 1), repmat(diff(G.cells.facePos)*nqc, nDof,1));
    WC = sparse(ii, jj, wc);
    
    wf = repmat(wf(:), nDof, 1);
    [bf, bc] = boundaryFaces(G);
    ncbf = sum((1:G.cells.num)' == bc',2);
    [ii, jj] = blockDiagIndex(ones(G.cells.num*nDof, 1), repmat((diff(G.cells.facePos) - ncbf)*nqf, nDof,1));
    WF = sparse(ii, jj, wf);
    
    % Accumulation term----------------------------------------------------
    
    sW_psi  = [];
    sW0_psi = [];
    if numel(pvMult) == 1
        pvMult = repmat(pvMult, G.cells.num,1);
    end
    if numel(pvMult0) == 1
        pvMult0 = repmat(pvMult0, G.cells.num,1);
    end
    
    for dofNo = 1:nDof
        sW_psi  = [sW_psi ; pvMult(cellNoc) .*bW(cellNoc) .*rock.poro(cellNoc).*sW(xc,cellNoc) .*psi{dofNo}(xc)];
        sW0_psi = [sW0_psi; pvMult0(cellNoc).*bW0(cellNoc).*rock.poro(cellNoc).*sW0(xc,cellNoc).*psi{dofNo}(xc)];
    end
    
    acc = WC*(sW_psi - sW0_psi)/dt;
    
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
    
    
    Gc = @(x,c) bW(c).*fW(x,c).*mobO(x,c);
    
    ig = [];
    for dofNo = 1:nDof
        ig = [ig;   fW(xc, cellNoc).*sum(vTc(cellNoc,:).*grad_psi{dofNo}(xc),2)  ...
                  + fW(xc, cellNoc).*mobO(xc, cellNoc).*sum((Gwc(cellNoc,:) - Goc(cellNoc,:)).*grad_psi{dofNo}(xc),2)];
    end
    flux1 = -WC*ig;
    
    faces = reshape(repmat(faces, nqf, 1), [], 1);
    
    
    upCells_v = G.faces.neighbors(:,2);
    
    intf = find(op.internalConn);
    upCells_v(intf(upcW)) = op.N(:,1);
    
%     upCells_v = rldecode(upCells_v, (diff(G.cells.facePos) - ncbf)*nqf);
    upCells_v = upCells_v(faces);
    
    
    upCells_G = upCells_v;
    
    ig = [];
    
    for dofNo = 1:nDof
        ig = [ig; (fW(xf, upCells_v).*vT(faces) + fW(xf, upCells_G).*mobO(xf, upCells_G).*(Gw(faces) - Go(faces))).*psi{dofNo}(xf)];
    end
    flux2 = WF*ig;
    
    flux = flux1 + flux2;
    
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

        [xwc, wwc, nqwc, ii, jj, cellNowc] = makeCellIntegrator(G, wc, model.degree, 'tri');
        xwc = (xwc - G.cells.centroids(cellNowc,:))./G.cells.diameters(cellNowc);
        
        wwc = repmat(wwc(:), nDof, 1);
        ncf = diff(G.cells.facePos);
        [ii, jj] = blockDiagIndex(ones(numel(wc)*nDof, 1), repmat(ncf(wc)*nqwc, nDof,1));
        WWC = sparse(ii, jj, wwc);
        
        ig = [];
        for dofNo = 1:nDof
            ig = [ig; bW(cellNowc).*wflux(cellNowc).*(fW(xwc, cellNowc).*(~isInj(cellNowc)) ...
                                                  + compPerf(cellNowc,1).*(isInj(cellNowc))).*psi{dofNo}(xwc)];
        end
        
        prod = WWC*ig;
        
        ind = mcolon((wc-1)*nDof + 1, wc*nDof);
        
        water(ind) = water(ind) - prod;
%         water = water - prod;
%         bWqW = bW(wc).*f_w_w.*wflux;
%         bOqO = bO(wc).*f_o_w.*wflux;

        % Store well fluxes
%         wflux_O = double(bOqO);
%         wflux_W = double(bWqW);
% 
%         for i = 1:numel(W)
%             perfind = perf2well == i;
%             state.wellSol(i).qOs = sum(wflux_O(perfind));
%             state.wellSol(i).qWs = sum(wflux_W(perfind));
%         end

    end

% % Get total flux from state
% flux = sum(state.flux, 2);
% vT = flux(model.operators.internalConn);
% 
% % Stored upstream indices
% [flag_v, flag_g] = getSaturationUpwind(model.upwindType, state, {Gw, Go}, vT, s.T, {mobW, mobO}, s.faceUpstr);
% flag = flag_v;
% 
% upcw  = flag(:, 1);
% upco  = flag(:, 2);
% 
% upcw_g = flag_g(:, 1);
% upco_g = flag_g(:, 2);
% 
% mobOf = s.faceUpstr(upco, mobO);
% mobWf = s.faceUpstr(upcw, mobW);
% 
% totMob = (mobOf + mobWf);
%     
% mobWf_G = s.faceUpstr(upcw_g, mobW);
% mobOf_G = s.faceUpstr(upco_g, mobO);
% mobTf_G = mobWf_G + mobOf_G;
% f_g = mobWf_G.*mobOf_G./mobTf_G;
% if opt.solveForWater
%     f_w = mobWf./totMob;
%     bWvW   = s.faceUpstr(upcw, bW).*f_w.*vT + s.faceUpstr(upcw_g, bO).*f_g.*s.T.*(Gw - Go);
% 
%     wat = (s.pv/dt).*(pvMult.*bW.*sW - pvMult0.*f.bW(p0).*sW0) + s.Div(bWvW);
%     if ~isempty(W)
%         wat(wc) = wat(wc) - bWqW;
%     end
% 
%     eqs = {wat};
%     names = {'water'};    
% else
%     f_o = mobOf./totMob;
%     bOvO   = s.faceUpstr(upco, bO).*f_o.*vT + s.faceUpstr(upco_g, bO).*f_g.*s.T.*(Go - Gw);
% 
%     oil = (s.pv/dt).*( pvMult.*bO.*(1-sW) - pvMult0.*f.bO(p0).*(1-sW0) ) + s.Div(bOvO);
%     if ~isempty(W)
%         oil(wc) = oil(wc) - bOqO;
%     end
%     eqs = {oil};
%     names = {'oil'};
% end
% types = {'cell'};
% rho = {rhoW, rhoO};
% mob = {mobW, mobO};
% sat = {sW, sO};
% [eqs, ~, src] = addBoundaryConditionsAndSources(model, eqs, names, types, state, ...
%                                      {pFlow, pFlow}, sat, mob, rho, ...
%                                      {}, {}, ...
%                                      drivingForces);

eqs = {water};
names = {'water'};
types = {'cell'};

% pv = reshape(repmat(op.pv', nDof, 1), [], 1);
% if ~model.useCNVConvergence
%     eqs{1} = eqs{1}.*(dt./pv);
% end

problem = LinearizedProblem(eqs, types, names, primaryVars, state, dt);
end

% function fun = getSatFromDof(sdof, psi, G)
%     
%     fun = @(x) x(:,1)*0;
%     nDof = size(psi,2);
%     for dofNo = 1:size(psi,2)
%         ix = (1:nDof:G.cells.num*nDof) + dofNo - 1;
%         fun = @(x) fun(x) + sdof(ix).*psi{dofNo}(x);
%     end
%     
% end