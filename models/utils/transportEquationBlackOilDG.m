function [problem, state] = transportEquationBlackOilDG(state0, state, model, dt, drivingForces, varargin)

    opt = struct('Verbose'      , mrstVerbose, ...
                 'reverseMode'  , false      , ...
                 'scaling'      , []         , ...
                 'resOnly'      , false      , ...
                 'solveForWater', false      , ...
                 'solveForOil'  , true       , ...
                 'solveForGas'  , true       , ...
                 'iteration'    , -1         , ...
                 'stepOptions'  , []         ); % Compatibility only
    opt      = merge_options(opt, varargin{:});
    
    % Frequently used properties
    op       = model.operators;
    fluid    = model.fluid;
    rock     = model.rock;
    G        = model.G;
    disc     = model.disc;
    flux2Vel = disc.velocityInterp.faceFlux2cellVelocity;
    W        = drivingForces.W;
    disgas   = model.disgas;
    vapoil   = model.vapoil;
    
    % We may solve for both oil and water simultaneously
    solveAllPhases = opt.solveForWater && opt.solveForOil && opt.solveForGas;
    
    % Prepare state for simulation-----------------------------------------
    if opt.iteration == 1 && ~opt.resOnly 
        if model.tryMaxDegree
            % If we are at the first iteration, we try to solve using
            % maximum degree in all cells
            state.degree(~G.cells.ghost) = disc.degree;
        end
        if ~isempty(W)
            state.degree(vertcat(W.cells)) = 0;
        end
        % For cells that previously had less than nDof unknowns, we must
        % map old dofs to new
        state = disc.mapDofs(state, state0, 's');
        state = disc.mapDofs(state, state0, 'rs');
        state = disc.mapDofs(state, state0, 'rv');
        
    end
    % Update discretizaiton information. This is carried by the state
    % variable, and holds the number of dofs per cell + dof position in
    % state.sdof
    state0 = disc.updateDofPos(state0);
    state  = disc.updateDofPos(state);
    %----------------------------------------------------------------------
    
    % Properties from current and previous timestep------------------------
    [p , sWdof , sOdof , sGdof , rSdof , rVdof , wellSol] = model.getProps(state , ...
                  'pressure', 'swdof', 'sodof', 'sgdof', 'rsdof', 'rvdof', 'wellsol');
    [p0, sWdof0, sOdof0, sGdof0, rSdof0, rVdof0         ] = model.getProps(state0, ...
                  'pressure', 'swdof', 'sodof', 'sgdof', 'rsdof', 'rvdof'           );
    % If timestep has been split relative to pressure, linearly interpolate
    % in pressure.
    if isfield(state, 'timestep')
        dt_frac = dt/state.timestep;
        p       = p.*dt_frac + p0.*(1-dt_frac);
    end
    %----------------------------------------------------------------------
    
    % Initialization of independent variables -----------------------------
    
    [sW , sO , sG ] = model.getProps(state , 'water', 'oil', 'gas');
    [sW0, sO0, sG0] = model.getProps(state0, 'water', 'oil', 'gas');
    st  = model.getCellStatusVO(state , sO ,  sW , sG );
    st0 = model.getCellStatusVO(state0, sO0,  sW0, sG0);
    
    if disgas || vapoil
        st  = cellfun(@(st) expandSingleValue(st, G), st , 'unif', false);
        st0 = cellfun(@(st) expandSingleValue(st, G), st0, 'unif', false);
        st  = cellfun(@(st) rldecode(st, state.nDof, 1), st, 'unif', false);
        st0 = cellfun(@(st) rldecode(st, state.nDof, 1), st0, 'unif', false);
        xDof = st{1}.*rSdof + st{2}.*rVdof + st{3}.*sGdof;
        gVar = 'xDof';
    else
        gVar = 'sGdof';
    end
    
    assert(~opt.reverseMode, 'Backwards solver not supported for splitting');
    if solveAllPhases
        if ~opt.resOnly
            if disgas || vapoil
                [sWdof, sOdof, xDof] = model.AutoDiffBackend.initVariablesAD(sWdof, sOdof, xDof);
            else
                [sWdof, sOdof, sGdof] = model.AutoDiffBackend.initVariablesAD(sWdof, sOdof, sGdof);
            end
        end
        primaryVars = {'sWdof', 'sOdof', gVar};
        sTdof = sOdof + sWdof + sGdof;
    else
        if ~opt.resOnly
            if disgas || vapoil
                [sWdof, xDof] = model.AutoDiffBackend.initVariablesAD(sWdof, xDof);
            else
                [sWdof, sGdof] = model.AutoDiffBackend.initVariablesAD(sWdof, sGdof);
            end
        end
        primaryVars = {'sWdof', gVar};
        sOdof     = -(sWdof + sGdof);
        ix        = disc.getDofIx(state, 1, Inf);
        sOdof(ix) = 1 + sOdof(ix);
        sTdof     = zeros(size(double(sWdof)));
        sTdof(ix) = 1;
    end
    
    if disgas || vapoil
        sGdof = st{2}.*(-sWdof);
        ix    = disc.getDofIx(state, 1, Inf);
        sGdof(ix) = sGdof(ix) + 1;
        sGdof = sGdof + st{3}.*xDof;
        if disgas
            rsSat = rldecode(fluid.rsSat(p), state.nDof, 1);
            rSdof = (~st{1}).*rsSat + st{1}.*xDof;
        end
        if vapoil
            rvSat = rldecode(fluid.rvSat(p), state.nDof, 1);
            rVdof = (~st{2}).*rvSat + st{2}.*xDof;
        end
    end
    
    %----------------------------------------------------------------------

    % Pressure and saturation dependent properties-------------------------
    % Get multipliers
    [pvMult, transMult, mobMult, pvMult0] = getMultipliers(model.fluid, p, p0);
    pvMult  = expandSingleValue(pvMult , G);
    pvMult0 = expandSingleValue(pvMult0, G);
    mobMult = expandSingleValue(mobMult, G);
    T       = op.T.*transMult;
    tm = ones(G.faces.num,1);
    tm(op.internalConn) = transMult;
    T_all   = model.operators.T_all.*tm;
    
    % Phase properties
    gdz = model.getGravityGradient();
    [b, mu, rho, mob] = getDerivedPropertyFunctionsBO(model, p, mobMult, st);
    bW   = b{1};   bO   = b{2};   bG   = b{3};
    rhoW = rho{1}; rhoO = rho{2}; rhoG = rho{3};
    mobW = mob{1}; mobO = mob{2}; mobG = mob{3};

    [sW, sO, sG] = disc.getCellMean(state, sWdof, sOdof, sGdof);
    
    % Gravity flux
    
    xf = G.faces.centroids(disc.internalConn);
    cL = disc.N(:,1);
    cR = disc.N(:,2);
    [sWL, sOL, sGL, rsL, rvL] = disc.evaluateDGVariable(xf, cL, state, sWdof, sOdof, sGdof, rSdof, rVdof);
    [sWR, sOR, sGR, rsR, rvR] = disc.evaluateDGVariable(xf, cR, state, sWdof, sOdof, sGdof, rSdof, rVdof); 
    
    gW = (rhoW(cL, sWL) + rhoW(cR, sWR))/2.*gdz;
    if isfield(fluid, 'pcOW')
        gW = gW - op.Grad(fluid.pcOW(sW));
    end
    gO = (rhoO(cL, rsL) + rhoO(cR, rsR))/2.*gdz;
    gG = (rhoG(cL, sGL, rvL) + rhoG(cR, sGR, rvR))/2.*gdz;
    if isfield(fluid, 'pcOG')
        gG = gG + op.Grad(fluid.pcOG(sG));
    end
    P = sparse(find(op.internalConn), 1:nnz(op.internalConn), 1, G.faces.num, nnz(op.internalConn));
    gW = P*gW; gO = P*gO; gG = P*gG;
    % Add gravity flux where we have BCs to get correct cell values
    bc = drivingForces.bc;
    if ~isempty(bc)
        error('bcs not implemented');
        BCcells = sum(G.faces.neighbors(bc.face,:), 2);
        dz = G.cells.centroids(BCcells, :) - G.faces.centroids(bc.face,:);
        g = model.getGravityVector();
        rhoWBC = rhoW(BCcells);
        rhoOBC = rhoO(BCcells);
        rhoGBC = rhoO(BCcells);
        gW(bc.face) = rhoWBC.*(dz*g');
        gO(bc.face) = rhoOBC.*(dz*g');
    end    
    TgW  = T_all.*gW;          % Gravity volumetric flux, water
    TgO  = T_all.*gO;          % Gravity volumetric flux, oil
    TgG  = T_all.*gG;          % Gravity volumetric flux, gas
    TgWc = flux2Vel(TgW);      % Map to cell velocity
    TgOc = flux2Vel(TgO);
    TgGc = flux2Vel(TgG);
    TgW  = TgW./G.faces.areas; % Convert to gravity flux
    TgO  = TgO./G.faces.areas;
    TgG  = TgG./G.faces.areas;
    
    % Viscous flux
    flux  = sum(state.flux,2);   % Total volumetric flux
    vT    = flux./G.faces.areas; % Total flux 
    vTc   = flux2Vel(flux);      % Map face fluxes to cell velocities
%     
%     [qb_c, qb_f, r_c, r_fg] = computeSequentialFluxesDG(disc, model, state, T, T_all, ...
%         {gW, gO, gG}, {mobW, mobO, mobG}, {bW, bO, bG}, {sWdof, sOdof, sGdof, sTdof}, {0*rSdof, rSdof, rVdof});
    
    % Well contributions---------------------------------------------------
    if ~isempty(W)
        % Total well flux, composition and mappings
        perf2well = getPerforationToWellMapping(W);
        wc        = vertcat(W.cells);
        wflux     = zeros(G.cells.num,1);
        wflux(wc) = sum(vertcat(wellSol.flux), 2)./G.cells.volumes(wc);
        isInj     = wflux > 0;
        compWell  = vertcat(W.compi);
        compPerf  = zeros(G.cells.num, 3);
        compPerf(wc,:) = compWell(perf2well,:);
        
        % Saturations at cubature points
        [~, xcw, wcNo] = disc.getCubature(wc, 'volume');
        [sW_w, sO_w, sG_w, sT_w, rS_w, rV_w] = disc.evaluateDGVariable(xcw, wcNo, state, sWdof, sOdof, sGdof, sTdof, rSdof, rVdof);
        
        mobWW = mobW(wcNo, sW_w, sT_w);
        mobOW = mobO(wcNo, sO_w, sT_w, rS_w);
        mobGW = mobG(wcNo, sG_w, sT_w, rV_w);
        mobTW = mobWW + mobOW + mobGW;
        
        fWW = ~isInj(wcNo).*sT_w.*mobWW./mobTW + isInj(wcNo).*compPerf(wcNo,1);
        fOW = ~isInj(wcNo).*sT_w.*mobOW./mobTW + isInj(wcNo).*compPerf(wcNo,2);
        fGW = ~isInj(wcNo).*sT_w.*mobGW./mobTW + isInj(wcNo).*compPerf(wcNo,3);
        
        fWW = ~isInj(wcNo).*fWW + isInj(wcNo).*compPerf(wcNo,1);
        fOW = ~isInj(wcNo).*fOW + isInj(wcNo).*compPerf(wcNo,2);
        fGW = ~isInj(wcNo).*fGW + isInj(wcNo).*compPerf(wcNo,3);
        
        bWW = bW(wcNo, sW_w);
        bOW = bO(wcNo, rS_w);
        bGW = bG(wcNo, sG_w, rV_w);
        
        qWW = wflux(wcNo).*sT_w.*fWW;
        qOW = wflux(wcNo).*sT_w.*fOW;
        qGW = wflux(wcNo).*sT_w.*fGW;
        
        % Water well contributions
        integrandW = @(psi, gradPsi) qWW.*bWW.*psi;
        srcW_w = disc.cellInt(integrandW, wc, state, sWdof);
        
        % Oil well contributions
        integrandO = @(psi, gradPsi) (qOW.*bOW + rV_w.*qGW.*bGW).*psi;
        srcO_w = disc.cellInt(integrandO, wc, state, sWdof);
        
        % Oil well contributions
        integrandG = @(psi, gradPsi) (qGW.*bGW + rS_w.*qOW.*bOW).*psi;
        srcG_w = disc.cellInt(integrandG, wc, state, sWdof);
                                 
        % Store well fluxes
        ix     = disc.getDofIx(state, 1, wc);
        wfluxW = double(srcW_w(ix));
        wfluxO = double(srcO_w(ix));
        wfluxG = double(srcG_w(ix));
        for wNo = 1:numel(W)
            perfind = perf2well == wNo;
            state.wellSol(wNo).qWs = sum(wfluxW(perfind));
            state.wellSol(wNo).qOs = sum(wfluxO(perfind));
            state.wellSol(wNo).qGs = sum(wfluxG(perfind));
        end

    end
    %----------------------------------------------------------------------

    % Evaluate saturation at cubature points-------------------------------
    % Cell cubature points
    [~, xc, c] = disc.getCubature((1:G.cells.num)', 'volume');
    [sWc , sOc , sGc , sTc , rvc , rsc] = disc.evaluateDGVariable(xc, c, state , sWdof , sOdof , sGdof , sTdof , rVdof , rSdof);
    [sWc0, sOc0, sGc0, rvc0, rsc0]      = disc.evaluateDGVariable(xc, c, state0, sWdof0, sOdof0, sGdof0, rVdof0, rSdof0);
    % Face cubature points
    [~, xf, ~, f] = disc.getCubature((1:G.cells.num)', 'surface');
    % Upstream cells
    [~, ~, cfV, cfG] = disc.getSaturationUpwind(f, xf, T, flux, state, ...
                            {gW, gO, gG}, {mobW, mobO, mobG}, {sWdof, sOdof, sGdof}, {0*rSdof, rSdof, rVdof});
    % Water saturation
    [sWfV , sTWfV] = disc.evaluateDGVariable(xf, cfV(:,1), state, sWdof, sTdof);
    [sWfG , sTWfG] = disc.evaluateDGVariable(xf, cfG(:,1), state, sWdof, sTdof);
    % Oil saturation
    [sOfV , sTOfV, rsfV] = disc.evaluateDGVariable(xf, cfV(:,2), state, sOdof, sTdof, rSdof);
    [sOfG , sTOfG, rsfG] = disc.evaluateDGVariable(xf, cfG(:,2), state, sOdof, sTdof, rSdof);
    % Gas saturation
    [sGfV , sTGfV, rvfV] = disc.evaluateDGVariable(xf, cfV(:,3), state, sGdof, sTdof, rVdof);
    [sGfG , sTGfG, rvfG] = disc.evaluateDGVariable(xf, cfG(:,3), state, sGdof, sTdof, rVdof);
    
    %----------------------------------------------------------------------
    
    mobWc = mobW(c, sWc, sTc);
    mobOc = mobO(c, sOc, sTc, rvc);
    mobGc = mobG(c, sGc, sTc, rsc);
    mobTc = mobWc + mobOc + mobGc;
    
    mobWfV = mobW(cfV(:,1), sWfV, sTWfV);
    mobOfV = mobO(cfV(:,2), sOfV, sTOfV);
    mobGfV = mobG(cfV(:,3), sGfV, sTGfV);
    mobTfV = mobWfV + mobOfV + mobGfV;
    
    mobWfG = mobW(cfG(:,1), sWfG, sTWfG);
    mobOfG = mobO(cfG(:,2), sOfG, sTOfG);
    mobGfG = mobG(cfG(:,3), sGfG, sTGfG);
    mobTfG = mobWfG + mobOfG + mobGfG;
    
    [eqs, names, types] = deal(cell(1,2 + solveAllPhases));
    [types{:}] = deal('cell');
    
    eqNo = 1;
    % Water equation-------------------------------------------------------
    if opt.solveForWater
        % Cell values
        fWc  = sTc.*mobWc./mobTc;
        bWc  = bW(c, sWc);
        bWc0 = bW(c, sWc0);
         
        fWfV = mobWfV./mobTfV;
        fWfG = mobWfG./mobTfG;
        bWfV = bW(cfV(:,1), sWfV, sTWfV);
        bWfG = bW(cfV(:,1), sWfG, sTWfG);
        
        % Accumulation term
        acc = @(psi) (pvMult(c) .*rock.poro(c).*bWc .*sWc - ...
                      pvMult0(c).*rock.poro(c).*bWc0.*sWc0).*psi/dt;
        % Convection term
        conv = @(gradPsi) bWc.*fWc.*(disc.dot(vTc(c,:) ...
                            + mobOc.*(TgWc(c,:) - TgOc(c,:)) ...
                            + mobGc.*(TgWc(c,:) - TgGc(c,:)),gradPsi));
        integrand = @(psi, gradPsi) acc(psi) - conv(gradPsi);
        % Integrate integrand*psi{dofNo} over all cells for dofNo = 1:nDof
        cellIntegralW = disc.cellInt(integrand, [], state, sWdof);
        % Flux term
        integrand = @(psi) ...
            (sTWfV.*bWfV.*fWfV.*vT(f) + bWfG.*fWfG.*( ...
                            mobOfG.*(TgW(f) - TgO(f)) ...
                          + mobGfG.*(TgW(f) - TgG(f)))).*psi;
        % Integrate integrand*psi{dofNo} over all cells surfaces for dofNo = 1:nDof
        faceIntegralW = disc.faceFluxInt(integrand, [], state, sWdof);
        % Sum integrals
        water = cellIntegralW + faceIntegralW;
        % Add well contributions
        if ~isempty(W)
            ix = disc.getDofIx(state, Inf, wc);
            water(ix) = water(ix) - srcW_w(ix);
        end
        eqs{eqNo}   = water;
        names{eqNo} = 'water';
        eqNo        = eqNo + 1;
    end
    %----------------------------------------------------------------------
    
    % Oil equation---------------------------------------------------------
    if opt.solveForOil
        % Cell values
        
        fOc  = sTc.*mobOc./mobTc;
        bOc  = bO(c, rvc);
        bOc0 = bO(c, rvc0);
         
        fOfV = mobOfV./mobTfV;
        fOfG = mobOfG./mobTfG;
        bOfV = bO(cfV(:,2), sWfV, sTWfV);
        bOfG = bO(cfV(:,2), sWfG, sTWfG);
        
        % Accumulation term
        acc = @(psi) (pvMult(c) .*rock.poro(c).*bOc .*sOc - ...
                      pvMult0(c).*rock.poro(c).*bOc0.*sOc0).*psi/dt;
        % Convection term
        conv = @(gradPsi) bOc.*fOc.*(disc.dot(vTc(c,:) ...
                            + mobWc.*(TgOc(c,:) - TgWc(c,:)) ...
                            + mobGc.*(TgOc(c,:) - TgGc(c,:)),gradPsi));
        integrand = @(psi, gradPsi) acc(psi) - conv(gradPsi);
        % Integrate integrand*psi{dofNo} over all cells for dofNo = 1:nDof
        cellIntegralO = disc.cellInt(integrand, [], state, sOdof);
        % Flux term
        integrand = @(psi) ...
            (sTOfV.*bOfV.*fOfV.*vT(f) + bOfG.*fOfG.*( ...
                            mobWfG.*(TgO(f) - TgW(f)) ...
                          + mobGfG.*(TgO(f) - TgG(f)))).*psi;
        % Integrate integrand*psi{dofNo} over all cells surfaces for dofNo = 1:nDof
        faceIntegralO = disc.faceFluxInt(integrand, [], state, sOdof);
        % Sum integrals
        oil = cellIntegralO + faceIntegralO;
        % Add well contributions
        if ~isempty(W)
            ix = disc.getDofIx(state, Inf, wc);
            oil(ix) = oil(ix) - srcO_w(ix);
        end
        eqs{eqNo}   = oil;
        names{eqNo} = 'oil';
        eqNo        = eqNo + 1;
    end
    %----------------------------------------------------------------------

    % Gas equation---------------------------------------------------------
    if opt.solveForGas
        % Cell values
        
        fGc  = sTc.*mobGc./mobTc;
        bGc  = bG(c, sGc);
        bGc0 = bG(c, sGc0);
         
        fGfV = mobGfV./mobTfV;
        fGfG = mobGfG./mobTfG;
        bGfV = bG(cfV(:,3), sGfV, sTGfV);
        bGfG = bG(cfV(:,3), sGfG, sTGfG);
        
        % Accumulation term
        acc = @(psi) (pvMult(c) .*rock.poro(c).*bGc .*sGc - ...
                      pvMult0(c).*rock.poro(c).*bGc0.*sGc0).*psi/dt;
        % Convection term
        conv = @(gradPsi) bGc.*fGc.*(disc.dot(vTc(c,:) ...
                            + mobWc.*(TgGc(c,:) - TgWc(c,:)) ...
                            + mobOc.*(TgGc(c,:) - TgOc(c,:)),gradPsi));
        integrand = @(psi, gradPsi) acc(psi) - conv(gradPsi);
        % Integrate integrand*psi{dofNo} over all cells for dofNo = 1:nDof
        cellIntegralG = disc.cellInt(integrand, [], state, sOdof);
        % Flux term
        integrand = @(psi) ...
            (sTGfV.*bGfV.*fGfV.*vT(f) + bGfG.*fGfG.*( ...
                            mobWfG.*(TgG(f) - TgW(f)) ...
                          + mobOfG.*(TgG(f) - TgO(f)))).*psi;
        % Integrate integrand*psi{dofNo} over all cells surfaces for dofNo = 1:nDof
        faceIntegralG = disc.faceFluxInt(integrand, [], state, sOdof);
        % Sum integrals
        gas = cellIntegralG + faceIntegralG;
        % Add well contributions
        if ~isempty(W)
            ix = disc.getDofIx(state, Inf, wc);
            gas(ix) = gas(ix) - srcG_w(ix);
        end
        eqs{eqNo}   = gas;
        names{eqNo} = 'gas';
    end
    %----------------------------------------------------------------------
    
    % Add BCs--------------------------------------------------------------
    if ~isempty(bc)
        % Boundary faces
        faces = bc.face;
        % Saturation outside boundary
        sL  = bc.sat;
        % Face cubature foordinates
        [~, x, ~, fBC] = disc.getCubature(faces, 'face');
        cBC = sum(G.faces.neighbors(fBC,:),2);
%         [xR, ~, ~] = disc.transformCoords(x, cBC);
        % Mapping from BC faces to global faces
        globFace2BCface        = nan(G.faces.num,1);
        globFace2BCface(faces) = 1:numel(faces);        
        locFaceNo = globFace2BCface(fBC);
        % Determine injectng boundaries
        sgn = 1 - 2*(G.faces.neighbors(faces, 1) == 0);
        isInj = vT(faces) > 0 & sgn < 0;
        % Upstream saturation
        [sW_r, sOR, sTR]  = disc.evaluateDGVariable(x, cBC, state, sWdof, sOdof, sTdof);
        sWBC = sL(locFaceNo,1).*isInj(locFaceNo) + sW_r.*(~isInj(locFaceNo));
        sOBC = sL(locFaceNo,2).*isInj(locFaceNo) + sOR.*(~isInj(locFaceNo));
        sTBC = sum(sL(locFaceNo,:),2).*isInj(locFaceNo) + sTR.*(~isInj(locFaceNo));
        % Frational flow functions
        fWBC = fW(sWBC, sOBC, sTBC, cBC, cBC);
        fOBC = fO(sWBC, sOBC, sTBC, cBC, cBC);
        if opt.solveForWater
            % Add water bc flux to water equation
            faceIntegrand = @(psi) (bW(cBC).*fWBC.*vT(fBC) ...
                      + bW(cBC).*fWBC.*mobO(sOBC,sTBC,cBC).*(TgW(fBC) - TgO(fBC))).*psi;
            fluxWBC = disc.faceFluxIntBC(faceIntegrand, bc, state, sWdof);
            water   = water + fluxWBC;
        end
        if opt.solveForOil
            % Add oil bc flux to oil equation
            faceIntegrand = @(psi) (bO(cBC).*fOBC.*vT(fBC) ...
                      + bO(cBC).*fOBC.*mobW(sWBC,sTBC,cBC).*(TgO(fBC) - TgW(fBC))).*psi;
            fluxOBC = disc.faceFluxIntBC(faceIntegrand, bc, state, sOdof);
            oil     = oil + fluxOBC;
        end
    end
    %----------------------------------------------------------------------
    
    % Add sources----------------------------------------------------------
    src = drivingForces.src;
    if ~isempty(src) 
        % Cubature
        [~, ~, cSRC] = disc.getCubature(src.cell, 'volume');
        % Mapping from source cells to global cells
        globCell2SRCcell = nan(G.cells.num,1);
        globCell2SRCcell(src.cell) = 1:numel(src.cell);
        cSRCloc = globCell2SRCcell(cSRC);
        % Total rate and saturaion at cubature points
        qT   = src.rate(cSRCloc)./G.cells.volumes(cSRC);
        sSRC = src.sat(cSRCloc,:);
        if opt.solveForWater
            % Add water source to water equation
            srcIntegrand = @(psi, gradPsi) bW(cSRC).*qT.*sSRC(:,1).*psi;
            srcW  = disc.cellInt(srcIntegrand, src.cell, state, sWdof);
            water = water - srcW;
        end
        if opt.solveForOil
            % Add oil source to oil equation
            srcIntegrand = @(psi, gradPsi) bO(cSRC).*qT.*sSRC(:,2).*psi;
            srcO = disc.cellInt(srcIntegrand, src.cell, state, sOdof);
            oil  = oil - srcO;
        end
    end
    %----------------------------------------------------------------------
    
    % Make Linearized problem----------------------------------------------
    % Scale equations
    if ~model.useCNVConvergence
        pv = rldecode(op.pv, state.nDof, 1);
        for eqNo = 1:numel(eqs)
            eqs{eqNo} = eqs{eqNo}.*(dt./pv);
        end
    end
    % Extra state output
    if model.extraStateOutput
        state     = model.storeDensity(state, rhoW, rhoO, []);
        state.cfl = dt.*sum(abs(vTc)./G.cells.dx,2);
    end
    % Linearize
    problem = LinearizedProblem(eqs, types, names, primaryVars, state, dt);
    %----------------------------------------------------------------------
    
    % Extract subproblem if we are solving subproblem----------------------
    if any(strcmpi(G.type, 'subgrid'))
        ix = disc.getDofIx(state, Inf, ~G.cells.ghost);
        
        for eqNo = 1:numel(problem.equations)
            eq = problem.equations{eqNo};
            if isa(eq, 'ADI')
                eq.val = eq.val(ix);
                eq.jac = cellfun(@(j) j(ix,ix), eq.jac, 'unif', false);
            else
                eq = eq(ix);
            end
            problem.equations{eqNo} = eq;
        end
    end
    %----------------------------------------------------------------------
    
end

% Expang single scalar values to one per cell------------------------------
function v = expandSingleValue(v,G)
    if numel(double(v)) == 1
        v = v*ones(G.cells.num,1);
    end
end
%--------------------------------------------------------------------------