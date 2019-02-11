function [problem, state] = transportEquationOilWaterDG(state0, state, model, dt, drivingForces, varargin)

    opt = struct('Verbose'      , mrstVerbose, ...
                 'reverseMode'  , false      , ...
                 'scaling'      , []         , ...
                 'resOnly'      , false      , ...
                 'history'      , []         , ...
                 'solveForWater', true       , ...
                 'solveForOil'  , false      , ...
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
    
    % We may solve for both oil and water simultaneously
    solveAllPhases = opt.solveForWater && opt.solveForOil;
    
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
        state = disc.mapDofs(state, state0);
        
    end
    % Update discretizaiton information. This is carried by the state
    % variable, and holds the number of dofs per cell + dof position in
    % state.sdof
    state0 = disc.updateDofPos(state0);
    state  = disc.updateDofPos(state);
    %----------------------------------------------------------------------
    
    % Properties from current and previous timestep------------------------
    [p , sWdof , sOdof , wellSol] = model.getProps(state , ...
                  'pressure', 'swdof', 'sodof', 'wellsol');
    [p0, sWdof0, sOdof0,        ] = model.getProps(state0, ...
                  'pressure', 'swdof', 'sodof'           );
    % If timestep has been split relative to pressure, linearly interpolate
    % in pressure.
    if isfield(state, 'timestep')
        dt_frac = dt/state.timestep;
        p       = p.*dt_frac + p0.*(1-dt_frac);
    end
    %----------------------------------------------------------------------
    
    % Initialization of independent variables -----------------------------
    assert(~opt.reverseMode, 'Backwards solver not supported for splitting');
    if solveAllPhases
        if ~opt.resOnly
            [sWdof, sOdof] = model.AutoDiffBackend.initVariablesAD(sWdof, sOdof);
        end
        primaryVars = {'sWdof', 'sOdof'};
        sTdof = sOdof + sWdof;
    else
        if ~opt.resOnly
            sWdof = model.AutoDiffBackend.initVariablesAD(sWdof);
        end
        primaryVars = {'sWdof'};
        sOdof     = -sWdof;
        ix        = disc.getDofIx(state, 1, Inf);
        sOdof(ix) = 1 + sOdof(ix);
        sTdof     = zeros(size(double(sWdof)));
        sTdof(ix) = 1;
    end
    %----------------------------------------------------------------------

    % Pressure and saturation dependent properties-------------------------
    % Get multipliers
    [pvMult, transMult, mobMult, pvMult0] = getMultipliers(model.fluid, p, p0);
    pvMult  = expandSingleValue(pvMult , G);
    pvMult0 = expandSingleValue(pvMult0, G);
    mobMult = expandSingleValue(mobMult, G);
    T       = op.T.*transMult;
    T_all   = model.operators.T_all;
    
    % Phase properties
    gdz = model.getGravityGradient();
    [b, mu, rho, mob] = getDerivedPropertyFunctionsBO(model, p, mobMult, []);
    bW = b{1}; bO = b{2};
    muW = mu{1}; muO = mu{2};
    rhoW = rho{1}; rhoO = rho{2};
    mobW = mob{1}; mobO = mob{2};
    
%     [vW, bW, mobW, rhoW, pW, upcW, dpW, muW] ...
%                              = getPropsWater_DG(model, p, T, gdz, mobMult);
%     [vO, bO, mobO, rhoO, pO, upcO, dpO, muO] ...
%                              = getPropsOil_DG(model, p, T, gdz, mobMult  );
    bW0 = fluid.bW(p0);
    bO0 = fluid.bO(p0);
    
                           
    % Fractional flow functions
    fW = @(sW, sO, sT, cW, cO) mobW(cW, sW, sT)./(mobW(cW, sW, sT) + mobO(cO, sO, sT));
    fO = @(sW, sO, sT, cW, cO) mobO(cO, sO, sT)./(mobW(cW, sW, sT) + mobO(cO, sO, sT));
    
    xf = G.faces.centroids(disc.internalConn);
    cL = disc.N(:,1);
    cR = disc.N(:,2);
    sWL = disc.evaluateDGVariable(xf, cL, state, sWdof);
    sWR = disc.evaluateDGVariable(xf, cR, state, sWdof);
    
    gW = (rhoW(cL, sWL) + rhoW(cR, sWR))/2.*gdz;
    if isfield(fluid, 'pcOW')
        gW = gW - op.Grad(fluid.pcOW(sW));
    end
    gO = (rhoO(cL) + rhoO(cR))/2.*gdz;
    P = sparse(find(op.internalConn), 1:nnz(op.internalConn), 1, G.faces.num, nnz(op.internalConn));
    gW = P*gW;
    gO = P*gO;
    
    % Gravity flux
%     gp = op.Grad(p);
%     [gW, gO] = deal(zeros(G.faces.num, 1));
%     gW(op.internalConn) = gp - dpW;
%     gO(op.internalConn) = gp - dpO;
    % Add gravity flux where we have BCs to get correct cell values
    bc = drivingForces.bc;
    if ~isempty(bc)
        BCcells = sum(G.faces.neighbors(bc.face,:), 2);
        dz = G.cells.centroids(BCcells, :) - G.faces.centroids(bc.face,:);
        g = model.getGravityVector();
        rhoWBC = rhoW(BCcells);
        rhoOBC = rhoO(BCcells);
        gW(bc.face) = rhoWBC.*(dz*g');
        gO(bc.face) = rhoOBC.*(dz*g');
    end    
    TgW  = T_all.*gW;          % Gravity volumetric flux, water
    TgO  = T_all.*gO;          % Gravity volumetric flux, oil
    TgWc = flux2Vel(TgW);      % Map to cell velocity
    TgOc = flux2Vel(TgO);
    TgW  = TgW./G.faces.areas; % Convert to gravity flux
    TgO  = TgO./G.faces.areas;
    
    % Viscous flux
    flux  = sum(state.flux,2);   % Total volumetric flux
    vT    = flux./G.faces.areas; % Total flux 
    vTc   = flux2Vel(flux);      % Map face fluxes to cell velocities
    %----------------------------------------------------------------------
    
    % Well contributions---------------------------------------------------
    if ~isempty(W)
        % Total well flux, composition and mappings
        perf2well = getPerforationToWellMapping(W);
        wc        = vertcat(W.cells);
        wflux     = zeros(G.cells.num,1);
        wflux(wc) = sum(vertcat(wellSol.flux), 2)./G.cells.volumes(wc);
        isInj     = wflux > 0;
        compWell  = vertcat(W.compi);
        compPerf  = zeros(G.cells.num, 2);
        compPerf(wc,:) = compWell(perf2well,:);
        
        % Saturations at cubature points
        [~, xcw, wcNo] = disc.getCubature(wc, 'volume');
        [sWW, sOW, sTW] = disc.evaluateDGVariable(xcw, wcNo, state, sWdof, sOdof, sTdof);
        
        % Water well contributions
        integrand = @(psi, gradPsi) bW(wcNo, sWW).*wflux(wcNo)...
            .*(sTW.*fW(sWW, sOW, sTW, wcNo, wcNo).*(~isInj(wcNo)) + compPerf(wcNo,1).*isInj(wcNo)).*psi;
        srcWW = disc.cellInt(integrand, wc, state, sWdof);
        
        % Oil well contributions
        integrand = @(psi, gradPsi) bO(wcNo).*wflux(wcNo)...
            .*(sTW.*fO(sWW, sOW, sTW, wcNo, wcNo).*(~isInj(wcNo)) + compPerf(wcNo,2).*isInj(wcNo)).*psi;
        srcOW = disc.cellInt(integrand, wc, state, sOdof);
                
        % Store well fluxes
        ix     = disc.getDofIx(state, 1, wc);
        wfluxW = double(srcWW(ix));
        wfluxO = double(srcOW(ix));
        for wNo = 1:numel(W)
            perfind = perf2well == wNo;
            state.wellSol(wNo).qWs = sum(wfluxW(perfind));
            state.wellSol(wNo).qOs = sum(wfluxO(perfind));
        end

    end
    %----------------------------------------------------------------------

    % Evaluate saturation at cubature points-------------------------------
    % Cell cubature points
    [~, xc, c] = disc.getCubature((1:G.cells.num)', 'volume');
    [sWc , sOc , sTc] = disc.evaluateDGVariable(xc, c, state , sWdof , sOdof , sTdof);
    [sWc0, sOc0]      = disc.evaluateDGVariable(xc, c, state0, sWdof0, sOdof0);
    % Face cubature points
    [~, xf, ~, f] = disc.getCubature((1:G.cells.num)', 'surface');
    % Upstream cells
    [~, ~, cfV, cfG] = disc.getSaturationUpwind(f, xf, T, flux, state, ...
                            {gW, gO}, {mobW, mobO}, {sWdof, sOdof}, {});
    % Water saturation
    [sWfV , sTWfV] = disc.evaluateDGVariable(xf, cfV(:,1), state, sWdof, sTdof);
    [sWfG , sTWfG] = disc.evaluateDGVariable(xf, cfG(:,1), state, sWdof, sTdof);
    % Oil saturation
    [sOfV , sTOfV] = disc.evaluateDGVariable(xf, cfV(:,2), state, sOdof, sTdof);
    [sOfG , sTOfG] = disc.evaluateDGVariable(xf, cfG(:,2), state, sOdof, sTdof);
    %----------------------------------------------------------------------
    
    % Water equation-------------------------------------------------------
    if opt.solveForWater
        % Cell values
        fWc   = fW(sWc,sOc,sTc,c,c);
        mobOc = mobO(c, sOc,sTc);
        % Face values
        fWfV   = fW(sWfV, sOfV, sTWfV, cfV(:,1), cfV(:,2));
        fWfG   = fW(sWfG, sOfG, sTWfG, cfG(:,1), cfG(:,2));
        mobOfG = mobO(cfG(:,2), sOfG, sTOfG);
        % Accumulation term
        acc = @(psi) (pvMult(c) .*rock.poro(c).*bW(c) .*sWc - ...
                      pvMult0(c).*rock.poro(c).*bW0(c).*sWc0).*psi/dt;
        % Convection term
        conv = @(gradPsi) bW(c).*fWc.*(disc.dot(vTc(c,:),gradPsi) ...
                        + mobOc.*disc.dot(TgWc(c,:) - TgOc(c,:),gradPsi));
        integrand = @(psi, gradPsi) acc(psi) - conv(gradPsi);
        % Integrate integrand*psi{dofNo} over all cells for dofNo = 1:nDof
        cellIntegralW = disc.cellInt(integrand, [], state, sWdof);
        % Flux term
        integrand = @(psi) ...
            (sTWfV.*bW(cfV(:,1)).*fWfV.*vT(f) ...
                  + bW(cfG(:,1)).*fWfG.*mobOfG.*(TgW(f) - TgO(f))).*psi;
        % Integrate integrand*psi{dofNo} over all cells surfaces for dofNo = 1:nDof
        faceIntegralW = disc.faceFluxInt(integrand, [], state, sWdof);
        % Sum integrals
        water = cellIntegralW + faceIntegralW;
        % Add well contributions
        if ~isempty(W)
            ix = disc.getDofIx(state, Inf, wc);
            water(ix) = water(ix) - srcWW(ix);
        end
        
    end
    %----------------------------------------------------------------------
    
    % Oil equation---------------------------------------------------------
    if opt.solveForOil
         % Cell values
        fOc   = fO(sWc,sOc,sTc,c,c);
        mobWc = mobW(c, sWc,sTc);
        % Face values
        fOfV   = fO(sWfV, sOfV, sTOfV, cfV(:,1), cfV(:,2));
        fOfG   = fO(sWfG, sOfG, sTOfG, cfG(:,1), cfG(:,2));
        mobWfG = mobW(cfG(:,1), sWfG, sTWfG);
        % Accumulation term
        acc = @(psi) (pvMult(c) .*rock.poro(c).*bO(c) .*sOc - ...
                      pvMult0(c).*rock.poro(c).*bO0(c).*sOc0).*psi/dt;
        % Convection term
        conv = @(gradPsi) bO(c).*fOc.*(disc.dot(vTc(c,:),gradPsi) ...
                        + mobWc.*disc.dot(TgOc(c,:) - TgWc(c,:),gradPsi));
        integrand = @(psi, gradPsi) acc(psi) - conv(gradPsi);
        % Integrate integrand*psi{dofNo} over all cells for dofNo = 1:nDof
        cellIntegralO = disc.cellInt(integrand, [], state, sOdof);
        % Flux term
        integrand = @(psi) ...
            (sTOfV.*bO(cfV(:,2)).*fOfV.*vT(f) ...
                  + bO(cfG(:,2)).*fOfG.*mobWfG.*(TgO(f) - TgW(f))).*psi;
        % Integrate integrand*psi{dofNo} over all cells surfaces for dofNo = 1:nDof
        faceIntegralO = disc.faceFluxInt(integrand, [], state, sOdof);
        % Sum integrals
        oil = cellIntegralO + faceIntegralO;
        % Add well contributions
        if ~isempty(W)
            ix = disc.getDofIx(state, Inf, wc);
            oil(ix) = oil(ix) - srcOW(ix);
        end
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
        [sWR, sOR, sTR]  = disc.evaluateDGVariable(x, cBC, state, sWdof, sOdof, sTdof);
        sWBC = sL(locFaceNo,1).*isInj(locFaceNo) + sWR.*(~isInj(locFaceNo));
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
    % Define equations, names and types
    if solveAllPhases
        eqs = {water, oil};
        names = {'water', 'oil'};    
        types = {'cell', 'cell'};
    elseif opt.solveForWater
        eqs   = {water  };
        names = {'water'};
        types = {'cell' };
    elseif opt.solveForOil
        eqs   = {oil   };
        names = {'oil' };
        types = {'cell'};
    end
    % Scale equations
    if ~model.useCNVConvergence
        pv = rldecode(op.pv, state.nDof, 1);
        for eqNo = 1:(opt.solveForOil+opt.solveForWater)
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
    if any(strcmpi(G.type, 'subgrid'));
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