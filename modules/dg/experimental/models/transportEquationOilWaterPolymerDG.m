function [problem, state] = transportEquationOilWaterPolymerDG(state0, state, model, dt, drivingForces, varargin)
%Undocumented Utility Function

%{
Copyright 2009-2023 SINTEF Digital, Mathematics & Cybernetics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}

    opt = struct('Verbose'        , mrstVerbose, ...
                 'reverseMode'    , false      , ...
                 'scaling'        , []         , ...
                 'resOnly'        , false      , ...
                 'history'        , []         , ...
                 'solveForWater'  , false      , ...
                 'solveForOil'    , true       , ...
                 'solveForPolymer', true       , ...
                 'iteration'      , -1         , ...
                 'stepOptions'    , []         ); % Compatibility only
    opt      = merge_options(opt, varargin{:});
    
    % Frequently used properties
    op       = model.operators;
    fluid    = model.fluid;
    rock     = model.rock;
    G        = model.G;
    disc     = model.disc;
    W        = drivingForces.W;
    psi      = disc.basis.psi;
    gradPsi  = disc.basis.grad_psi;
    
    % We may solve for both oil and water simultaneously
    solveAllPhases = opt.solveForWater && opt.solveForOil;
    
    % Prepare state for simulation-----------------------------------------
    if opt.iteration == 1 && ~opt.resOnly 
        if model.tryMaxDegree
            % If we are at the first iteration, we try to solve using
            % maximum degree in all cells
            state.degree(~G.cells.ghost) = disc.degree;
        end
        if ~isempty(disc.degree0)
            state.degree = disc.degree0;
        end
        if ~isempty(W)
            state.degree(vertcat(W.cells)) = 0;
        end
        % For cells that previously had less than nDof unknowns, we must
        % map old dofs to new 
        state = disc.mapDofs(state, state0, 's');
        state = disc.mapDofs(state, state0, 'c');
        
    end
    % Update discretizaiton information. This is carried by the state
    % variable, and holds the number of dofs per cell + dof position in
    % state.sdof
    state0        = disc.updateDofPos(state0);
    [state, disc] = disc.updateDofPos(state);
    %----------------------------------------------------------------------
    
    % Properties from current and previous timestep------------------------
    [p , sWdof , sOdof , cDof , cMax , wellSol] = model.getProps(state , ...
                  'pressure', 'swdof', 'sodof', 'cdof', 'polymermax', 'wellsol');
    [p0, sWdof0, sOdof0, cDof0, cMax0         ] = model.getProps(state0, ...
                  'pressure', 'swdof', 'sodof', 'cdof', 'polymermax'           );
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
            [sWdof, sOdof, cDof] = model.AutoDiffBackend.initVariablesAD(sWdof, sOdof, cDof);
        end
        primaryVars = {'sWdof', 'sOdof', 'cDof'};
        sTdof = sOdof + sWdof;
    else
        if ~opt.resOnly
            [sWdof, cDof] = model.AutoDiffBackend.initVariablesAD(sWdof, cDof);
        end
        primaryVars = {'sWdof', 'cDof'};
        sOdof     = -sWdof;
        ix        = disc.getDofIx(state, 1, Inf);
        sOdof(ix) = 1 + sOdof(ix);
        sTdof     = zeros(size(double(sWdof)));
        sTdof(ix) = 1;
    end
    disc.sample = sWdof;
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
    [bW, bO, muW, muO, rhoW, rhoO, mobW, mobO] = deal(b{:}, mu{:}, rho{:}, mob{:});
    [b0, mu0, rho0, mob0] = getDerivedPropertyFunctionsBO(model, p0, mobMult, []);
    [bW0, bO0, muW0, muO0, rhoW0, rhoO0, mobW0, mobO0] = deal(b0{:}, mu0{:}, rho0{:}, mob0{:});
    
    [muWMult, a, cbar] = getMobilityMultipliers(model, cMax);
    mobW = @(e, s, r, c) mobW(e, s)./muWMult(e, c);
    
    x_f = G.faces.centroids(disc.internalConn);
    c_l = disc.N(:,1);
    c_r = disc.N(:,2);
    sW_l = disc.evaluateDGVariable(x_f, c_l, state, sWdof);
    sW_r = disc.evaluateDGVariable(x_f, c_r, state, sWdof);
    sW = disc.getCellMean(state, sWdof);
    % Gravity flux
    gW = (rhoW(c_l, sW_l) + rhoW(c_r, sW_r))/2.*gdz;
    if isfield(fluid, 'pcOW')
        gW = gW - op.Grad(fluid.pcOW(sW));
    end
    gO = (rhoO(c_l) + rhoO(c_r))/2.*gdz;
    P = sparse(find(op.internalConn), 1:nnz(op.internalConn), 1, G.faces.num, nnz(op.internalConn));
    gW = P*gW;
    gO = P*gO;
    
    % Add gravity flux where we have BCs to get correct cell values
    bc = drivingForces.bc;
    if ~isempty(bc)
        [bv_bc, g] = computeBoundaryFluxesDG(model, state, bc, T_all, ...
            {gW, gO}, {mobW, mobO}, {bW, bO}, {rhoW, rhoO}, {sWdof, sOdof});
        [bWvW_bc, bOvO_bc] = deal(bv_bc{:});
        [gW, gO] = deal(g{:});
    end
    
    [bv_c, bv_f, ~, bWvP_c, bWvP_f] ...
        = computeSequentialFluxesDG(disc, model, state, T, T_all, ...
        {gW, gO}, {mobW, mobO}, {bW, bO}, {sWdof, sOdof, sTdof}, {}, cDof);
    [bWvW_c, bOvO_c] = deal(bv_c{:});
    [bWvW_f, bOvO_f] = deal(bv_f{:});
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
        [~, x_w, e_w] = disc.getCubature(wc, 'volume');
        [sW_w, sO_w, sT_w, c_w] = disc.evaluateDGVariable(x_w, e_w, state, sWdof, sOdof, sTdof, cDof);
        s_w = {sW_w./sT_w, sO_w./sT_w};
        
        mobW_w = mobW(e_w, s_w, 0, c_w);
        mobO_w = mobO(e_w, s_w);
        mobT_w = mobW_w + mobO_w;
        
        fW_w = ~isInj(e_w).*sT_w.*mobW_w./mobT_w + isInj(e_w).*compPerf(e_w,1);
        fO_w = ~isInj(e_w).*sT_w.*mobO_w./mobT_w + isInj(e_w).*compPerf(e_w,2);
        
        bWqW_w = bW(e_w, sW_w).*wflux(e_w).*sT_w.*fW_w;
        bOqO_w = bO(e_w, sO_w).*wflux(e_w).*sT_w.*fO_w;
        
        % Polymer well equations
        c_winj = getWellPolymerDG(G, W, perf2well, e_w);
        c_w = c_w.*(~isInj(e_w)) + c_winj.*(isInj(e_w));
        bWqP_w  = c_w.*bWqW_w;
        
        % Water well contributions
        srcW_w = disc.inner(bWqW_w, psi, 'dV', wc);
        % Oil well contributions
        srcO_w = disc.inner(bOqO_w, psi, 'dV', wc);
        % Polymer well contributions
        srcP_w = disc.inner(bWqP_w, psi, 'dV', wc);
        
        % Store well fluxes
        ix     = disc.getDofIx(state, 1, wc);
        wfluxW = double(srcW_w(ix));
        wfluxO = double(srcO_w(ix));
        wfluxP = double(srcP_w(ix));
        for wNo = 1:numel(W)
            perfind = perf2well == wNo;
            state.wellSol(wNo).qWs = sum(wfluxW(perfind));
            state.wellSol(wNo).qOs = sum(wfluxO(perfind));
            state.wellSol(wNo).qPs = sum(wfluxP(perfind));
        end

    end
    %----------------------------------------------------------------------

    % Evaluate saturation at cubature points-------------------------------
    % Cell cubature points
    [~, x_e, e] = disc.getCubature((1:G.cells.num)', 'volume');
    [sW_c , sO_c , sT_c] = disc.evaluateDGVariable(x_e, e, state , sWdof , sOdof , sTdof);
    [sW0_c, sO0_c]      = disc.evaluateDGVariable(x_e, e, state0, sWdof0, sOdof0);
    % B-factors at current timestep
    bW_c  = bW(e, sW_c);
    bO_c  = bO(e, sO0_c);
    % B-factors at previous timestep
    bW0_c = bW0(e, sW0_c);
    bO0_c = bO0(e, sO0_c);
    
    poro = rock.poro;
    %----------------------------------------------------------------------
    
    [eqs, names, types] = deal(cell(2 + solveAllPhases,1));
    [types{:}] = deal('cell');
    eqNo = 1;
    % Water equation-------------------------------------------------------
    if opt.solveForWater
        % Mass terms
        mW  = pvMult(e) .*poro(e).*bW_c .*sW_c;
        mW0 = pvMult0(e).*poro(e).*bW0_c.*sW0_c;
        % Assemble inner products
        water =   disc.inner((mW - mW0)/dt, psi    , 'dV') ...
                - disc.inner(bWvW_c       , gradPsi, 'dV') ...
                + disc.inner(bWvW_f       , psi    , 'dS');
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
        mO  = pvMult(e) .*poro(e).*(bO_c .*sO_c);
        mO0 = pvMult0(e).*poro(e).*(bO0_c.*sO0_c);
        % Assmeble inner products
        oil =   disc.inner((mO - mO0)/dt, psi    , 'dV') ...
              - disc.inner(bOvO_c       , gradPsi, 'dV') ...
              + disc.inner(bOvO_f       , psi    , 'dS');
        % Add well contributions
        if ~isempty(W)
            ix = disc.getDofIx(state, Inf, wc);
            oil(ix) = oil(ix) - srcO_w(ix);
        end
        eqs{eqNo}   = oil;
        names{eqNo} = 'oil';
        eqNo = eqNo + 1;
    end
    %----------------------------------------------------------------------

    % Polymer equation-----------------------------------------------------
    if opt.solveForPolymer
        % Cell values
        c_e    = disc.evaluateDGVariable(x_e, e, state, cDof);
        c0_e   = disc.evaluateDGVariable(x_e, e, state0, cDof0);
        ads_e  = effads(c_e, cMax(e), model);
        ads0_e = effads(c0_e, cMax0(e), model);
        % Mass terms
        mP  = (1-fluid.dps).*pvMult(e) .*poro(e).*bW_c .*sW_c .*c_e;
        mP0 = (1-fluid.dps).*pvMult0(e).*poro(e).*bW0_c.*sW0_c.*c0_e;
        % Adsorption terms
        ads  = fluid.rhoR.*(1-poro(e)).*ads_e;
        ads0 = fluid.rhoR.*(1-poro(e)).*ads0_e;
        % Assmeble inner products
        polymer =   disc.inner((mP  - mP0 )/dt, psi    , 'dV') ...
                  + disc.inner((ads - ads0)/dt, psi    , 'dV') ...
                  - disc.inner(bWvP_c         , gradPsi, 'dV') ...
                  + disc.inner(bWvP_f         , psi    , 'dS');
        % Add well contributions
        if ~isempty(W)
            ix = disc.getDofIx(state, Inf, wc);
            polymer(ix) = polymer(ix) - srcP_w(ix);
        end
        % Fix model weakness
        if isa(polymer, 'ADI')
            isPolymer = strcmpi(primaryVars, 'cDof');
            epsilon = 1.e-8;
            epsilon = sqrt(epsilon)*mean(abs(diag(polymer.jac{isPolymer})));
            bad     = abs(diag(polymer.jac{isPolymer})) < epsilon;
            polymer(bad) = cDof(bad);
        end
        bad = double(sW) == 0;
        if any(bad)
            ix          = disc.getDofIx(state, Inf, bad);
            polymer(ix) = cDof(ix);
        end
        eqs{eqNo}   = polymer;
        names{eqNo} = 'polymer';
    end
    %----------------------------------------------------------------------
    
    % Add BCs--------------------------------------------------------------
    if ~isempty(bc)
        if opt.solveForWater
            % Add water bc flux to water equation
            eqs{eqNo} = eqs{eqNo} + disc.inner(bWvW_bc, psi, 'dSbc', [], bc);
            eqNo = eqNo + 1;
        end
        if opt.solveForOil
            % Add oil bc flux to oil equation
            eqs{eqNo} = eqs{eqNo} + disc.inner(bOvO_bc, psi, 'dSbc', [], bc);
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
    isPolymer = strcmpi(names, 'polymer');
    eqs{isPolymer} = eqs{isPolymer}./fluid.cmax;
    if ~model.useCNVConvergence
        pv = rldecode(op.pv, state.nDof, 1);
        for eqNo = 1:(2 + solveAllPhases)
            eqs{eqNo} = eqs{eqNo}.*(dt./pv);
        end
    end
    % Extra state output
    if model.extraStateOutput
        state     = model.storeDensity(state, rhoW, rhoO, []);
        state.cfl = dt.*sum(abs(vTc)./G.cells.dx,2);
    end
    % Linearize
    state.cmax0 = cMax0;
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

    if 0
        eqsDG = eqs;
        save('dgEqs.mat', 'eqsDG');
    end
    
end

% Expang single scalar values to one per cell------------------------------
function v = expandSingleValue(v,G)
    if numel(double(v)) == 1
        v = v*ones(G.cells.num,1);
    end
end
%--------------------------------------------------------------------------

% Effective adsorption, depending of desorption or not
function y = effads(c, cmax, model)
   if model.fluid.adsInx == 2
      y = model.fluid.ads(max(c, cmax));
   else
      y = model.fluid.ads(c);
   end
end
%--------------------------------------------------------------------------

% Multipliers due to polymer-----------------------------------------------
function [muWMult, a, cbar] = getMobilityMultipliers(model, cMax)
    fluid = model.fluid;
    ads = @(e, c) effads(c, cMax(e), model);
    mixpar = fluid.mixPar;
    cbar   = @(c) c/fluid.cmax;
    a = fluid.muWMult(fluid.cmax).^(1-mixpar);
    b = @(c) 1./(1-cbar(c)+cbar(c)./a);
    % The viscosity multiplier only result from the polymer mixing.
    muWeffMult = @(c) b(c).*fluid.muWMult(c).^mixpar;
    permRed = @(e, c) 1 + ((fluid.rrf-1)./fluid.adsMax).*ads(e, c);
    muWMult = @(e, c) muWeffMult(c).*permRed(e, c);
end
%--------------------------------------------------------------------------
