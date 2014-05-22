function [problem, state] = equationsGasWater(state0, state, dt, G, drivingForces, s, f, t, varargin)

    opt = struct('Verbose'     , mrstVerbose , ...
                 'reverseMode' , false       , ...
                 'resOnly'     , false       , ...
                 'minerals'    , false       , ...
                 'iteration'   , -1          , ...
                 'stepOptions' , []); % compatibility only
    opt = merge_options(opt, varargin{:});
    
    
    W = drivingForces.Wells;
    assert(isempty(drivingForces.src)); 

    if isempty(W) && isfield(state, 'wellSol') && ~isempty(state.wellSol) 
        % no wells present - empty content of wellSol structure
        state.wellSol.qGs = [];
        state.wellSol.qWs = [];
        state.wellSol.bhp = [];
    end

    % ----------------------------- Current variables -----------------------------
    p   = state.pressure;
    sG  = state.s(:,2); % water is first column, gas second
    pBH = vertcat(state.wellSol.bhp);
    qWs = vertcat(state.wellSol.qWs);
    qGs = vertcat(state.wellSol.qGs);
    
    % ---------------------------- Previous variables ----------------------------
    p0  = state0.pressure;
    sG0 = state0.s(:,2);
        
    % ------------------ Initialization of independent variables ------------------
    
    if ~opt.resOnly
        if ~opt.reverseMode
            [p, sG, qWs, qGs, pBH] = initVariablesADI(p, sG, qWs, qGs, pBH);
        else
            [p0, sG0, ~, ~, ~] = initVariablesADI(p0, sw0, 0*qWs, 0*qGs, 0*pBH);
        end
        primaryVars = {'pressure', 'sG', 'qWs', 'qGs', 'bhp'};
    end
        
    % ----------------------------------------------------------------------------
    
    % Check for p-dependent tran mult:
    trMult = 1;
    if isfield(f, 'tranMultR'), trMult = f.tranMultR(p); end;
    
    % Check for p-dependent porv mult:
    pvMult = 1; pvMult0 = 1;
    if isfield(f, 'pvMultR')
        pvMult  = f.pvMultR(p);
        pvMult0 = f.pvMultR(p0);
    end
    transMult = 1;
    if isfield(f, 'transMult')
        transMult = f.transMult(p);
    end

    trans = s.T .* transMult;
    
    % Check for capillary pressure
    pcGW = 0;
    if isfield(f, 'pcGW')
        pcGW = f.pcGW(sG);
    end
    % ----------------------------------------------------------------------------
    
    [krW, krG] = f.relPerm(sG);
    
    % computing densities, mobilities and upstream indices
    [bW, mobW, fluxW] = compMFlux(p - pcGW, t, f.bW, f.muW, f.rhoWS, trMult, krW, s, G, trans);
    [bG, mobG, fluxG] = compMFlux(p,        t, f.bG, f.muG, f.rhoGS, trMult, krG, s, G, trans);

    bW0 = f.bW(p0 - pcGW, t);
    bG0 = f.bG(p0, t);
    
    % --------------------------- Continuity equations ---------------------------
    
    eqs{1} = (s.pv/dt) .* (pvMult .* bW .* (1-sG) - pvMult0 .* bW0 .* (1-sG0)) + s.div(fluxW);
    eqs{2} = (s.pv/dt) .* (pvMult .* bG .* sG     - pvMult0 .* bG0 .* sG0    ) + s.div(fluxG);
    names  = {'water', 'gas'};
    types  = {'cell' , 'cell'};
    
    % ---------------------------- Boundary conditions ----------------------------

    [bc_cells, b_fW, b_fG] = ...
        BCFluxes(G, s, p, t, f, drivingForces.bc, krW, krG, transMult, trMult);

    % @@ If a cell has more than one pressure boundary condition, only one
    % will be taken into account.  This should be fixed later.
    eqs{1}(bc_cells) = eqs{1}(bc_cells) + b_fW;
    eqs{2}(bc_cells) = eqs{2}(bc_cells) + b_fG;
    
    % ------------------------------ Well equations ------------------------------
    if ~isempty(W) && ~opt.reverseMode
    
        wc = vertcat(W.cells);
        rhos = [f.rhoWS, f.rhoGS];
        bw = {bW(wc), bG(wc)};
        pw = p(wc);
        mw = {mobW(wc), mobG(wc)};
        
        [eqs(3:5), cqs, state.wellSol] = ...
            getWellContributions(W, state.wellSol, pBH, {qWs, qGs}, pw, rhos, ...
                                 bw, {}, {}, mw, 'iteration', opt.iteration, 'model', 'WG'); 
        [wc, cqs] = checkForRepetitions(wc, cqs);
        eqs{1}(wc) = eqs{1}(wc) - cqs{1};
        eqs{2}(wc) = eqs{2}(wc) - cqs{2};

        names(3:5) = {'waterWells' , 'gasWells' , 'closureWells'};
        types(3:5) = {'perf'       , 'perf'     , 'well'};
    
    else % no wells, or in reverse mode
        eqs(3:5)   = {0*pBH   , 0*pBH   , 0*pBH}; % empty ADIs
        names(3:5) = {'empty' , 'empty' , 'empty'};
        types(3:5) = {'none'  , 'none'  , 'none'};
    end
    % ----------------------------------------------------------------------------
    
    problem = linearProblem(eqs, types, names, primaryVars, state);
    
end

% ============================= END MAIN FUNCTION =============================

function [cells, fluxW, fluxG] = BCFluxes(G, s, p, t, f, bc, krW, krG, transMult, trMult)
    
    if isempty(bc)
        % all boundary conditions are no-flow; nothing to do here.
        cells = []; fluxW = []; fluxG = []; return;
    end

    assert(all(strcmp(bc.type, 'pressure'))); % only supported type for now
    cells = sum(G.faces.neighbors(bc.face, :), 2);
    %assert(numel(unique(cells)) == numel(cells)); % multiple BC per cell not supported
    
    % prepare boundary-specific values
    bp   = p(cells); 
    bt   = t(cells);
    bdp  = bc.value - bp;
    bdz  = G.faces.centroids(bc.face, 3) - G.cells.centroids(cells,3);
    
    assert(isscalar(transMult)); % not implemented for face-wise values
    assert(isscalar(trMult));    % not implemented for face-wise values
    trans = s.T_all(bc.face) .* transMult; 
    krWf = trMult * krW(cells);
    krGf = trMult * krG(cells);
    
    g = norm(gravity); 
    
    [bw, rhoWf, mobW] = computeRhoMobBFace(bp, bt, f.bW, f.rhoWS, krWf, f.muW);
    [bg, rhoGf, mobG] = computeRhoMobBFace(bp, bt, f.bG, f.rhoGS, krGf, f.muG);
    
    dptermW = bdp - rhoWf .* g .* bdz;
    dptermG = bdp - rhoGf .* g .* bdz;
    
    % Adjust upstream-weighted mobilities to prevent gas from re-entering the domain
    ix = dptermG > 0;
    mobG(ix) = 0;
    mobW(ix) = trMult ./ f.muW(bp(ix), bt(ix));
 
    % compute fluxes
    fluxW = - bw .* mobW .* trans .* dptermW;
    fluxG = - bg .* mobG .* trans .* dptermG;
    
    %fprintf('%f, %f\n', max(abs(fluxW.val)), max(abs(fluxG.val)));
end

% ----------------------------------------------------------------------------

function [bb, brho, bmob] = computeRhoMobBFace(bp, bt, bfun, rhoS, krf, mufun)

    bb   = bfun(bp, bt);
    brho = bb .* rhoS;
    bmob = krf ./ mufun(bp, bt);

end

% ----------------------------------------------------------------------------
function [b, mob, flux] = compMFlux(p, t, bfun, mufun, rhoS, trMult, kr, s, G, trans)
    
    g   = norm(gravity);
    b   = bfun(p, t);
    mob = trMult .* kr ./ mufun(p, t);

    rhoFace = s.faceAvg(b*rhoS);
    
    dp   = s.grad(p) - g * rhoFace .* s.grad(G.cells.centroids(:,3));
    upc  = (double(dp)>=0);
    flux = s.faceUpstr(upc, b .* mob) .* trans .* dp;

end

% ----------------------------------------------------------------------------
function [wc, cqs] = checkForRepetitions(wc, cqs)
[c, ~, ic] = unique(wc, 'stable');
if numel(c) ~= numel(wc)
    A = sparse(ic, (1:numel(wc))', 1, numel(c), numel(wc));
    wc = c;
    for k=1:numel(cqs)
        cqs{k} = A*cqs{k};
    end
end
end
