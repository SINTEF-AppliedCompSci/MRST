function [problem, state] = transportEquationOilWaterDG(state0, state, model, dt, drivingForces, varargin)
%Undocumented Utility Function

%{
Copyright 2009-2020 SINTEF Digital, Mathematics & Cybernetics.

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
    state0        = disc.updateDofPos(state0);
    [state, disc] = disc.updateDofPos(state);
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
    
    [bv_c, bv_f] = computeSequentialFluxesDG(disc, model, state, T, T_all, ...
        {gW, gO}, {mobW, mobO}, {bW, bO}, {sWdof, sOdof, sTdof}, {});
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
        [~, x_w, c_w] = disc.getCubature(wc, 'volume');
        [sW_w, sO_w, sT_w] = disc.evaluateDGVariable(x_w, c_w, state, sWdof, sOdof, sTdof);
        s_w = {sW_w./sT_w, sO_w./sT_w};
        
        mobW_w = mobW(c_w, s_w);
        mobO_w = mobO(c_w, s_w);
        mobT_w = mobW_w + mobO_w;
        
        fW_w = ~isInj(c_w).*sT_w.*mobW_w./mobT_w + isInj(c_w).*compPerf(c_w,1);
        fO_w = ~isInj(c_w).*sT_w.*mobO_w./mobT_w + isInj(c_w).*compPerf(c_w,2);
        
        bWqW_w = bW(c_w, sW_w).*wflux(c_w).*sT_w.*fW_w;
        bOqO_w = bO(c_w, sO_w).*wflux(c_w).*sT_w.*fO_w;
        
        % Water well contributions
        srcW_w = disc.inner(bWqW_w, psi, 'dV', wc);
        % Oil well contributions
        srcO_w = disc.inner(bOqO_w, psi, 'dV', wc);
                
        % Store well fluxes
        ix     = disc.getDofIx(state, 1, wc);
        wfluxW = double(srcW_w(ix));
        wfluxO = double(srcO_w(ix));
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
    [sW_c , sO_c , sT_c] = disc.evaluateDGVariable(xc, c, state , sWdof , sOdof , sTdof);
    [sW0_c, sO0_c]      = disc.evaluateDGVariable(xc, c, state0, sWdof0, sOdof0);
    % B-factors at current timestep
    bW_c  = bW(c, sW_c);
    bO_c  = bO(c, sO0_c);
    % B-factors at previous timestep
    bW0_c = bW0(c, sW0_c);
    bO0_c = bO0(c, sO0_c);
    %----------------------------------------------------------------------
    
    [eqs, names, types] = deal(cell(opt.solveForWater + opt.solveForOil,1));
    [types{:}] = deal('cell');
    eqNo = 1;
    % Water equation-------------------------------------------------------
    if opt.solveForWater
        % Mass terms
        mW  = pvMult(c) .*rock.poro(c).*bW_c .*sW_c;
        mW0 = pvMult0(c).*rock.poro(c).*bW0_c.*sW0_c;
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
        mO  = pvMult(c) .*rock.poro(c).*(bO_c .*sO_c);
        mO0 = pvMult0(c).*rock.poro(c).*(bO0_c.*sO0_c);
        % Assmeble inner products
        oil =   disc.inner((mO - mO0)/dt, psi    , 'dV') ...
              - disc.inner(bOvO_c       , gradPsi, 'dV') ...
              + disc.inner(bOvO_f       , psi    , 'dS');
        % Add well contributions
        if ~isempty(W)
            ix = disc.getDofIx(state, Inf, wc);
            oil(ix) = oil(ix) - srcO_w(ix);
        end
        if disc.weakJumpLimiter
            oil = oil + disc.getJump([], state, sOdof);
        end
        eqs{eqNo}   = oil;
        names{eqNo} = 'oil';
    end
    %----------------------------------------------------------------------

    % Add BCs--------------------------------------------------------------
    eqNo = 1;
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
    if ~model.useCNVConvergence
        pv = rldecode(op.pv, state.nDof, 1);
        for eqNo = 1:(opt.solveForOil+opt.solveForWater)
            eqs{eqNo} = eqs{eqNo}.*(dt./pv);
        end
    end
    % Extra state output
    if model.extraStateOutput
%         state     = model.storeDensity(state, rhoW_c, rhoO_c, []);
%         state.cfl = dt.*sum(abs(vTc)./G.cells.dx,2);
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
