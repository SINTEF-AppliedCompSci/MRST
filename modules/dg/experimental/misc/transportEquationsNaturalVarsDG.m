function [problem, state] = transportEquationsNaturalVarsDG(state0, state, model, dt, drivingForces, varargin)
    opt = struct('Verbose'   , mrstVerbose,...
                'reverseMode', false      ,...
                'resOnly'    , false      ,...
                'iteration'  , -1         );
    opt = merge_options(opt, varargin{:});

    % Frequently used properties
    op       = model.operators;
    fluid    = model.fluid;
    rock     = model.rock;
    G        = model.G;
    disc     = model.disc;
    flux2Vel = disc.velocityInterp.faceFlux2cellVelocity;
    nDofMax  = disc.basis.nDof;

    % Prepare state for simulation-----------------------------------------
    if opt.iteration == 1 && ~opt.resOnly 
        if model.tryMaxDegree
            % If we are at the first iteration, we try to solve using
            % maximum degree in all cells
            state.degree(~G.cells.ghost) = disc.degree;
        end
        % For cells that previously had less than nDof unknowns, we must
        % map old dofs to new
        state = disc.mapDofs(state, state0, 's');
        state = disc.mapDofs(state, state0, 'x');
        state = disc.mapDofs(state, state0, 'y');
    end
    % Update discretizaiton information. This is carried by the state
    % variable, and holds the number of dofs per cell + dof position in
    % state.sdof
    state0 = disc.updateDofPos(state0);
    state  = disc.updateDofPos(state);
    
    oIx = 1 + model.water;
    gIx = 2 + model.water;
    
    ncomp = model.EOSModel.fluid.getNumberOfComponents();
    [x,y] = deal(zeros(size(state.x)));
    for cNo = 1:ncomp
        x(:,cNo) = model.disc.getCellMean(state.xdof(:,cNo), state);
        y(:,cNo) = model.disc.getCellMean(state.ydof(:,cNo), state);
    end
    dx = state.x./x;
    dx(isnan(dx)) = 1;
    state.xdof = state.xdof.*rldecode(dx, state.nDof, 1);
    for cNo = 1:size(state.x,2)
        ix = model.disc.getDofIx(state, 1, isinf(dx(:,cNo)));
        state.xdof(ix,cNo) = state.x(isinf(dx(:,cNo)), cNo);
        ix = model.disc.getDofIx(state, 2:nDofMax, isinf(dx(:,cNo)));
        state.xdof(ix,cNo) = 0;
    end
    dy = state.y./y;
    dy(isnan(dy)) = 1;
    state.ydof = state.ydof.*rldecode(dy, state.nDof, 1);
    for cNo = 1:size(state.y,2)
        ix = model.disc.getDofIx(state, 1, isinf(dy(:,cNo)));
        state.ydof(ix,cNo) = state.y(isinf(dy(:,cNo)), cNo);
        ix = model.disc.getDofIx(state, 2:nDofMax, isinf(dy(:,cNo)));
        state.ydof(ix,cNo) = 0;
    end
    s = zeros(size(state.s));
    for phNo = 1:size(state.s,2)
        s(:,phNo) = model.disc.getCellMean(state.sdof(:,phNo), state);
    end
    ds = state.s./s;
    ds(isnan(ds)) = 1;
    state.sdof = state.sdof.*rldecode(ds, state.nDof, 1);
    for phNo = 1:size(state.s,2)
        ix = model.disc.getDofIx(state, 1, isinf(ds(:,phNo)));
        state.sdof(ix,phNo) = state.s(isinf(ds(:,phNo)), phNo);
        ix = model.disc.getDofIx(state, 2:nDofMax, isinf(ds(:,phNo)));
        state.sdof(ix,phNo) = 0;
    end

%     fluid = model.fluid;
    compFluid = model.EOSModel.fluid;

    % Properties at current timestep
    [p, sWdof, sOdof, sGdof, xdof, ydof, temp, wellSol] = model.getProps(state, ...
        'pressure', 'swdof', 'sodof', 'sgdof', 'xdof', 'ydof', 'T', 'wellSol');
    assert(all(p>0), 'Pressure must be positive for compositional model');

    [p0, sWdof0, sOdof0, sGdof0, xdof0, ydof0, temp0, wellSol0] = model.getProps(state0, ...
        'pressure', 'swdof', 'sodof', 'sgdof', 'xdof', 'ydof', 'T', 'wellSol');

    if opt.iteration == 1
        state.change = zeros(G.cells.num,1);
    else
        state.change = state.change + (state.twoPhase ~= state0.twoPhase);
    end
    
    if any(state.change > 4) && 0
        pureLiquid = state0.pureLiquid;
        pureVapor  = state0.pureVapor;
        twoPhase   = state0.twoPhase;
    else
        [pureLiquid, pureVapor, twoPhase] = model.getFlag(state);
    end
    state.pureLiquid = pureLiquid;
    state.pureVapor  = pureVapor;
    state.twoPhase = twoPhase;
    
    pureLiquidIx = disc.getDofIx(state, Inf, pureLiquid);
    pureVaporIx  = disc.getDofIx(state, Inf, pureVapor);

    if 1
        % TODO: Do the same for dofs
        sO = state.s(:,1 + model.water);
        sG = state.s(:,2 + model.water);
        stol = 1e-8;
        pureWater = sO + sG < stol | sG == stol | sO == stol;
        sO(~pureVapor & pureWater) = stol;
        ix = disc.getDofIx(state, 1, ~pureVapor & pureWater);
        sOdof(ix) = stol;
        sG(~pureLiquid & pureWater) = stol;
        ix = disc.getDofIx(state, 1, ~pureLiquid & pureWater);
        sGdof(ix) = stol;

        sO0 = state0.s(:,1 + model.water);
        sG0 = state0.s(:,2 + model.water);
        [pureLiquid0, pureVapor0, twoPhase0] = model.getFlag(state0);
        pureWater0 = sO0 + sG0 < stol | sG == stol | sO == stol;
        sO0(~pureVapor0 & pureWater0) = stol;
        ix = disc.getDofIx(state, 1, ~pureVapor0 & pureWater0);
        sOdof0(ix) = stol;
        sG0(~pureLiquid0 & pureWater0) = stol;
        ix = disc.getDofIx(state, 1, ~pureLiquid0 & pureWater0);
        sGdof0(ix) = stol;
    end


    % if isfield(state, 'timestep') && opt.iteration == 1
    %     p = state.pressure_full;
    %     dt_frac = dt/state.timestep;
    %     state.pressure = p.*dt_frac + p0.*(1-dt_frac);
    % end

    zdof = state.componentsdof;
    ix = disc.getDofIx(state, Inf, ~twoPhase);
    xdof(ix, :) = zdof(ix, :);
    ydof(ix, :) = zdof(ix, :);
    if 1
        % TODO: Make function ensureMinimumFraction for dofs
        xdof = ensureMinimumFractionDG(disc, xdof, state);
        [ydof, z_tol] = ensureMinimumFractionDG(disc, ydof, state);
%         [y, z_tol] = ensureMinimumFraction(y);
    end
    xdof = expandMatrixToCell(xdof);
    ydof = expandMatrixToCell(ydof);

    [xnames, ynames, cnames] = deal(model.EOSModel.fluid.names);
    for i = 1:ncomp
        xnames{i} = ['v_', cnames{i}];
        ynames{i} = ['w_', cnames{i}];
    end
    twoPhaseIx = disc.getDofIx(state, Inf, twoPhase);
%     twoPhaseIx = find(twoPhase);

    wtmp = ones(sum(state.nDof(twoPhase)),1);
    wdof = cell(ncomp, 1);
    [wdof{:}] = deal(wtmp);

    nc = model.G.cells.num;
    for i = 1:(ncomp-1)
        wdof{i} = ydof{i}(twoPhaseIx);
    end
%     ix = disc.getDofIx(state, Inf, twoPhaseIx);
    sodof = sOdof(twoPhaseIx);

    X = sGdof;
    ix = disc.getDofIx(state, Inf, pureLiquid);
    X(ix,:) = sOdof(ix,:);

if model.water
    if ~opt.resOnly
        [X, xdof{1:ncomp-1}, sWdof, wdof{1:ncomp-1}, sodof] = model.AutoDiffBackend.initVariablesAD(...
         X, xdof{1:ncomp-1}, sWdof, wdof{1:ncomp-1}, sodof);
    end
    primaryVars = {'sGsO', xnames{1:end-1}, 'satw', ynames{1:end-1}, 'sato'};

else
    if ~opt.resOnly
        [X, xdof{1:ncomp-1}, wdof{1:ncomp-1}, sodof] = model.AutoDiffBackend.initVariablesAD(...
         X, xdof{1:ncomp-1}, wdof{1:ncomp-1}, sodof);    
    end
    primaryVars = {'sGsO', xnames{1:end-1}, ynames{1:end-1}, 'sato'};

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
    
sample = xdof{1};

sOdof = model.AutoDiffBackend.convertToAD(sOdof, sample);
sOdof(twoPhaseIx) = sodof;
sOdof(pureLiquidIx) = X(pureLiquidIx);
sGdof = model.AutoDiffBackend.convertToAD(sGdof, sample);
notPureLiquidIx = disc.getDofIx(state, Inf, ~pureLiquid);
sGdof(notPureLiquidIx) = X(notPureLiquidIx);

if model.water
    sTdof = sOdof + sGdof + sWdof;
    sTdof0 = sOdof0 + sGdof0 + sWdof0;
%     sWdof = sWdof./sTdof;
%     sWdof0 = sWdof0./sTdof0;
else
    sTdof = sOdof + sGdof;
    sTdof0 = sOdof0 + sGdof0;
end
% sT = max(sT, 1e-8);
% sT0 = max(sT0, 1e-8);

% sOdof = sOdof./sTdof;
% sGdof = sGdof./sTdof;

% sOdof0 = sOdof0./sTdof0;
% sGdof0 = sGdof0./sTdof0;


xdof{end} = zeros(size(double(xdof{1})));
ix = disc.getDofIx(state, 1);
xdof{end}(ix) = 1;


wdof{end} = zeros(size(double(wdof{1})));
ix = disc.getDofIx(state, 1, twoPhase);
P = sparse(twoPhaseIx, 1:numel(twoPhaseIx), 1, sum(state.nDof), numel(twoPhaseIx));
w = P*wdof{end};
w(ix) = 1;
wdof{end} = w(twoPhaseIx);
% wdof{end}(ix) = 1;

% xdof{end} = 1;
% wdof{end} = 1;
for cNo = 1:ncomp-1
    xdof{end} = xdof{end} - xdof{cNo};
    if any(twoPhase)
        wdof{end} = wdof{end}-wdof{cNo};
    end
end

for cNo = 1:ncomp
    ydof{cNo} = xdof{cNo};
    if any(twoPhase)
        ydof{cNo}(twoPhaseIx) = wdof{cNo};
    end
    xdof{cNo}(pureVaporIx)  = double(xdof{cNo}(pureVaporIx));
    ydof{cNo}(pureLiquidIx) = double(ydof{cNo}(pureLiquidIx));
end

cellJacMap = cell(numel(primaryVars), 1);
% toReorder = 1:nc;

% if isempty(twoPhaseIx) || opt.resOnly
%     reorder = [];
% else
%     % TODO: Implement for dg. Remember twoPhaseIx not the same as in
%     % original code
% %     tp = find(twoPhase);
%     nDof = sum(state.nDof);
%     n2ph = nnz(twoPhaseIx);
%     nVars = sum(sample.getNumVars());
%     reorder = 1:nVars;
%     start = twoPhaseIx + nDof;
%     stop = (nVars-n2ph+1):nVars;
%     
%     reorder(start) = stop;
%     reorder(stop) = start;
%     
%     offset = ncomp+model.water;
%     for i = 1:ncomp
%         cellJacMap{i + offset} = find(twoPhase);% twoPhaseIx;
%     end
% end
if isempty(twoPhaseIx) || opt.resOnly
    reorder = [];
else
    nDof = sum(state.nDof);
    n2ph = nnz(twoPhaseIx);
    nVars = sum(sample.getNumVars());
    reorder = 1:nVars;
    start = twoPhaseIx + nDof;
    stop = (nVars-n2ph+1):nVars;
    
    reorder(start) = stop;
    reorder(stop) = start;
    
    offset = ncomp+model.water;
    for i = 1:ncomp
        cellJacMap{i + offset} = find(twoPhase);
    end
    
%     cellno = rldecode(a, b);
    cellNo = rldecode((1:G.cells.num)', state.nDof, 1);
    
    for i = 2:offset
        cellJacMap{i} = cellNo;
    end
end


% Get cell averages of dg variables
[x, y] = deal(cell(ncomp,1));
z = zeros(model.G.cells.num, ncomp);
for cNo = 1:ncomp
    x{cNo}   = model.disc.getCellMean(xdof{cNo}, state);
    y{cNo}   = model.disc.getCellMean(ydof{cNo}, state);
    z(:,cNo) = model.disc.getCellMean(zdof(:,cNo), state);
end
% x = ensureMinimumFraction(x);
% y = ensureMinimumFraction(y);

sO = model.disc.getCellMean(sOdof, state);
sG = model.disc.getCellMean(sGdof, state);
sT = model.disc.getCellMean(sTdof, state);
sO = sO./sT;
sG = sG./sT;
if model.water
    sW = model.disc.getCellMean(sWdof, state);
    sW = sW./sT;
end

% Compute properties and fugacity
[xM,  yM,  rhoO,  rhoG,  muO,  muG, f_L, f_V, xM0, yM0, rhoO0, rhoG0] = ...
    model.getTimestepPropertiesEoS(state, state0, p, temp, x, y, z, sO, sG, cellJacMap);



% [pvMult, transMult, mobMult, pvMult0] = getMultipliers(model.fluid, p, p0);

if model.water
    sat = {sW, sO, sG};
else
    sat = {sO, sG};
end


% Compute transmissibility
T = transMult.*op.T;

rhoOf = op.faceAvg(sO.*rhoO)./max(op.faceAvg(sO), 1e-8);
rhoGf = op.faceAvg(sG.*rhoG)./max(op.faceAvg(sG), 1e-8);
mobO = @(sO,sT,c) mobMult(c).*fluid.krO(sO./sT)./muO(c);
mobG = @(sG,sT,c) mobMult(c).*fluid.krG(sG./sT)./muG(c);


% Gravity flux
gdz = model.getGravityGradient();
gO  = rhoOf.*gdz;
gG  = rhoGf.*gdz;
if isfield(fluid, 'pcOG')
    gG = gG - op.Grad(fluid.pcOG(sG));
end
P = sparse(find(op.internalConn), 1:nnz(op.internalConn), 1, G.faces.num, nnz(op.internalConn));
gO = P*gO;
gG = P*gG;

% g = 


% Viscous flux
flux  = sum(state.flux,2);   % Total volumetric flux
vT    = flux./G.faces.areas; % Total flux 
vTc   = flux2Vel(flux);      % Map face fluxes to cell velocities

if model.water
    pW    = p;
    pW0   = p0;
    bW    = fluid.bW(pW);
    rhoW  = bW.*fluid.rhoWS;
    rhoW0 = fluid.bW(pW0).*fluid.rhoWS;

    rhoWf = op.faceAvg(rhoW);
    muW   = fluid.muW(pW);
    mobW  = @(sW,sT,c) mobMult(c).*fluid.krW(sW./sT)./muW(c);

    gW    = rhoWf.*gdz;

    if isfield(fluid, 'pcOW')
        gW = gW + op.Grad(fluid.pcOW(sW));
    end
    gW = P*gW;
    sWt       = sW.*sT;
    g         = {gW, gO, gG};
    mob       = {mobW, mobO, mobG};
    rho       = {rhoW, rhoO, rhoG};
    pressures = {pW, p, p};
    sdof      = {sWdof, sOdof, sGdof};
else
    [rhoW, rhoW0, mobW, bW, sWt] = deal([]);
    g         = {gO, gG};
    mob       = {mobO, mobG};
    rho       = {rhoO, rhoG};
    sdof      = {sOdof, sGdof};
    pressures = {p, p};
end

% Add gravity flux where we have BCs to get correct cell values
bc = drivingForces.bc;
if ~isempty(bc)
    BCcells = sum(G.faces.neighbors(bc.face,:), 2);
    dz = G.cells.centroids(BCcells, :) - G.faces.centroids(bc.face,:);
    grav = model.getGravityVector();
    rhoWBC = rhoW(BCcells);
    rhoOBC = rhoO(BCcells);
    rhoGBC = rhoG(BCcells);
    gW(bc.face) = rhoWBC.*(dz*grav');
    gO(bc.face) = rhoOBC.*(dz*grav');
    gG(bc.face) = rhoGBC.*(dz*grav');
    g = {gW, gO, gG};
end

if model.water
TgW  = T_all.*gW;          % Gravity volumetric flux, water
TgWc = {disc.velocityInterp.D{1}*TgW, disc.velocityInterp.D{2}*TgW, disc.velocityInterp.D{3}*TgW};
TgWc = SpatialVector(TgWc{:});
TgW  = TgW./G.faces.areas; % Convert to gravity flux
end
TgO  = T_all.*gO;          % Gravity volumetric flux, oil
TgG  = T_all.*gG;          % Gravity volumetric flux, oil
% TgWc = flux2Vel(TgW);      % Map to cell velocity

TgOc = {disc.velocityInterp.D{1}*TgO, disc.velocityInterp.D{2}*TgO, disc.velocityInterp.D{3}*TgO};
TgOc = SpatialVector(TgOc{:});
TgGc = {disc.velocityInterp.D{1}*TgG, disc.velocityInterp.D{2}*TgG, disc.velocityInterp.D{3}*TgG};
TgGc = SpatialVector(TgGc{:});
% TgOc = flux2Vel(TgO);

TgO  = TgO./G.faces.areas;
TgG  = TgG./G.faces.areas;

    frac = fractionalFlowFunctionsDG(mob);
    
    % Well contributions---------------------------------------------------
    W = drivingForces.W;
    if ~isempty(W)
        % Total well flux, composition and mappings
        perf2well = getPerforationToWellMapping(W);
        wc        = vertcat(W.cells);
        [wflux, mflux] = deal(zeros(G.cells.num,1));
        wflux(wc) = sum(vertcat(wellSol.flux), 2)./G.cells.volumes(wc);
        mflux(wc) = sum(vertcat(wellSol.components), 2)./G.cells.volumes(wc);
        w_comp = vertcat(W.components);
        a = w_comp(perf2well, :).*repmat(compFluid.molarMass, numel(wc), 1);
        w_comp = zeros(G.cells.num,ncomp);
        w_comp(wc,:) = bsxfun(@rdivide, a, sum(a, 2));
        
        
        isInj     = wflux > 0;
        compWell  = vertcat(W.compi);
        compPerf  = zeros(G.cells.num, 2 + model.water);
        compPerf(wc,:) = compWell(perf2well,:);
        
        % Saturations at cubature points
        [xW, yW] = deal(cell(size(xdof)));
        [~, xicw, cNow] = disc.getCubature(wc, 'volume');
        [sWW, sOW, sGW, sTW, xW{:}, yW{:}] = disc.evaluateDGVariable(xicw, cNow, state, sWdof, sOdof, sGdof, sTdof, xdof{:}, ydof{:});
%         xW = ensureMinimumFraction(xW);
%         yW = ensureMinimumFraction(xW);
        
        if model.water
            fWW   = frac{1}({sWW,sOW,sGW}, sTW, [cNow, cNow, cNow]);
            fOW   = frac{2}({sWW,sOW,sGW}, sTW, [cNow, cNow, cNow]);
            fGW   = frac{3}({sWW,sOW,sGW}, sTW, [cNow, cNow, cNow]);
        else
            fOW   = frac{1}({sOW,sGW}, sTW, [cNow, cNow]);
            fGW   = frac{2}({sOW,sGW}, sTW, [cNow, cNow]);
        end

        xMW = model.EOSModel.getMassFraction(xW);
        xMW = cellfun(@(x) x.*sTW, xMW, 'unif', false);
        yMW = model.EOSModel.getMassFraction(yW);
        yMW = cellfun(@(y) y.*sTW, yMW, 'unif', false);
        
        srcW = cell(1, ncomp + model.water);
        if model.water
            % Water well contributions
            integrand = @(psi, gradPsi) rhoW(cNow).*wflux(cNow)...
                .*(sTW.*fWW.*(~isInj(cNow)) + compPerf(cNow,1).*isInj(cNow)).*psi;
            srcWW = disc.cellInt(integrand, wc, state, sWdof);
            srcW{ncomp + 1} = srcWW;
        end
        
        rOqO = rhoO(cNow).*fOW.*wflux(cNow).*(~isInj(cNow)) ...
               + compPerf(cNow, oIx).*mflux(cNow).*isInj(cNow);
        rGqG = rhoG(cNow).*fGW.*wflux(cNow) ...
               + compPerf(cNow, gIx).*mflux(cNow).*isInj(cNow);
%         rOqO(isInj(cNow)) = 
%         rGqG(isInj(cNow)) = compPerf(isInj(cNow), 2 + model.water).*mflux(isInj(cNow));
        
        for cNo = 1:ncomp
            integrand = @(psi, gradPsi) ...
            ( ...
                (rOqO + rGqG).*w_comp(cNow, cNo).*isInj(cNow)     ...
              + (rOqO.*xMW{cNo} + rGqG.*yMW{cNo}).*(~isInj(cNow)) ...
            ).*psi;
            srcW{cNo} = disc.cellInt(integrand, wc, state, sOdof);
        end
        
%         % Store well fluxes
%         ix     = disc.getDofIx(state, 1, wc);
%         wfluxW = double(srcWW(ix));
%         wfluxO = double(srcOW(ix));
%         for wNo = 1:numel(W)
%             perfind = perf2well == wNo;
%             state.wellSol(wNo).qWs = sum(wfluxW(perfind));
%             state.wellSol(wNo).qOs = sum(wfluxO(perfind));
%         end
        
        isProd = ~isInj;
        
        integrand =  @(psi, gradPsi) ...
            double(rOqO.*(sTW.*isProd(cNow) + ~isProd(cNow)))./fluid.rhoOS;
        qo = disc.cellInt(integrand, wc, state, state.sdof(:,1));
        
        integrand =  @(psi, gradPsi) ...
            double(rGqG.*(sTW.*isProd(cNow) + ~isProd(cNow)))./fluid.rhoGS;
        qg = disc.cellInt(integrand, wc, state, state.sdof(:,1));

        
%         qo(isProd) = qo(isProd);%.*sT(wc(isProd));
%         qg(isProd) = qg(isProd);%.*sT(wc(isProd));
        ix = disc.getDofIx(state, 1, wc);
        qo = qo(ix);
        qg = qg(ix);
        if model.water
            qw = double(srcWW(ix));
%              integrand =  @(psi, gradPsi) ...
%             double(rWqW.*(sTW.*isProd(cNow) + ~isProd(cNow)))./fluid.rhoWS;
%             qw = disc.cellInt(integrand, wc, state, state.sdof(:,1));
%             qw = double(rWqW)./fluid.rhoWS;
%             qw(isProd) = qw(isProd).*sT(wc(isProd));
        end
        
        compSrc = cell2mat(cellfun(@double, srcW, 'unif', false));
        compSrc = compSrc(ix,:);
        wf = wflux(wc).*G.cells.volumes(wc);
        for wNo = 1:numel(W)
            perfIx = perf2well == wNo;
            state.wellSol(wNo).components = compSrc(perfIx, :);
            state.wellSol(wNo).qGs        = sum(qg(perfIx));
            state.wellSol(wNo).qOs        = sum(qo(perfIx));
            if model.water
                state.wellSol(wNo).qWs = sum(qw(perfIx));
            end

            state.wellSol(wNo).qTr    = sum(wf(perfIx));
            state.wellSol(wNo).status = true;
            state.wellSol(wNo).type   = W(wNo).type;
            state.wellSol(wNo).val    = W(wNo).val;
        end
%         for i = 1:numel(W)
%             state.wellSol(i).components = compSrc(perf2well == i, :);
%             state.wellSol(i).qGs = sum(qg(perf2well == i));
%             state.wellSol(i).qOs = sum(qo(perf2well == i));
%             if model.water
%                 state.wellSol(i).qWs = sum(qw(perf2well == i));
%             end
% 
%             state.wellSol(i).qTr = sum(wflux(perf2well == i));
%             state.wellSol(i).status = true;
%             state.wellSol(i).type = W(i).type;
%             state.wellSol(i).val = W(i).val;
%         end
                
%         % Store well fluxes
%         ix     = disc.getDofIx(state, 1, wc);
%         wfluxW = double(srcWW(ix));
%         wfluxO = double(srcOW(ix));
%         for wNo = 1:numel(W)
%             perfind = perf2well == wNo;
%             state.wellSol(wNo).qWs = sum(wfluxW(perfind));
%             state.wellSol(wNo).qOs = sum(wfluxO(perfind));
%         end

    end

    % Evaluate saturation/compositions at cubature points------------------
    % Cell cubature points
    xdof0 = expandMatrixToCell(xdof0);
    ydof0 = expandMatrixToCell(ydof0);
    [xc, yc, xc0, yc0, xOfv, xGfv, yOfv, yGfv] = deal(cell(size(xdof)));
    
    [~, xic, c] = disc.getCubature((1:G.cells.num)', 'volume');
    [sWc , sOc , sGc , sTc , xc{:} , yc{:} ] = disc.evaluateDGVariable(xic, c, state , ...
                             sWdof , sOdof , sGdof , sTdof , xdof{:} , ydof{:});
%     xc = ensureMinimumFraction(xc);
%     yc = ensureMinimumFraction(yc);
    [sWc0, sOc0, sGc0, sTc0, xc0{:}, yc0{:}] = disc.evaluateDGVariable(xic, c, state0, ...
                             sWdof0, sOdof0, sGdof0, sTdof0, xdof0{:}, ydof0{:});
%     xc0 = ensureMinimumFraction(xc0);
%     yc0 = ensureMinimumFraction(yc0);
    % Face cubature pointsl
    [~, xif, ~, f] = disc.getCubature((1:G.cells.num)', 'surface');
    % Upstream cells
    [~, ~, cfv, cfg] = disc.getSaturationUpwind(f, xif, T, flux, g, mob, sdof, state);
    % Water saturation
    if model.water
        [sWfv , sTWfv] = disc.evaluateDGVariable(xif, cfv(:,1), state, sWdof, sTdof);
        [sWfg , sTWfg] = disc.evaluateDGVariable(xif, cfg(:,1), state, sWdof, sTdof);
    end
    
    % Oil saturation
    oind = 1 + model.water;
    [sOfv , sTOfv, xOfv{:}, yOfv{:}] = disc.evaluateDGVariable(xif, cfv(:,oind), state, sOdof, sTdof, xdof{:}, ydof{:});
%     xOfv = ensureMinimumFraction(xOfv);
%     yOfv = ensureMinimumFraction(yOfv);
    [sOfg , sTOfg] = disc.evaluateDGVariable(xif, cfg(:,oind), state, sOdof, sTdof);
    % Gas saturation
    gind = 2 + model.water;
    [sGfv , sTGfv, xGfv{:}, yGfv{:}] = disc.evaluateDGVariable(xif, cfv(:,gind), state, sGdof, sTdof, xdof{:}, ydof{:});
    [sGfg , sTGfg] = disc.evaluateDGVariable(xif, cfg(:,gind), state, sGdof, sTdof);
%     xGfv = ensureMinimumFraction(xGfv);
%     yGfv = ensureMinimumFraction(yGfv);
    
    %----------------------------------------------------------------------
    
    xMc   = model.EOSModel.getMassFraction(xc);
    xMc0  = model.EOSModel.getMassFraction(xc0);
    xMOfv = model.EOSModel.getMassFraction(xOfv);
    xMGfv = model.EOSModel.getMassFraction(xGfv);
    
    yMc   = model.EOSModel.getMassFraction(yc);
    yMc0  = model.EOSModel.getMassFraction(yc0);
    yMOfv = model.EOSModel.getMassFraction(yOfv);
    yMGfv = model.EOSModel.getMassFraction(yGfv);
    
%     for cNo = 1:ncomp
%        xMc{cNo} =  xMc{cNo}.*sTc;
%        yMc{cNo} =  yMc{cNo}.*sTc;
%        
%        xMc0{cNo} = xMc0{cNo}.*sTc0;
%        yMc0{cNo} = yMc0{cNo}.*sTc0;
%        
%        xMOfv{cNo} = xMOfv{cNo}.*sTOfv;
%        yMGfv{cNo} = yMGfv{cNo}.*sTGfv;
%     end
    
%     fWc   = frac{1}({sWc ,sOc ,sGc }, sTc  , [c,c,c]);
%     fWv   = frac{1}({sWfv,sOfv,sGfv}, sTWfv, cfv    );
%     fWg   = frac{1}({sWfg,sOfg,sGfg}, sTWfg, cfg    );

%     fOc   = frac{2}({sWc ,sOc ,sGc }, sTc  , [c,c,c]);
%     fOv   = frac{2}({sWfv,sOfv,sGfv}, sTOfv, cfv    );
%     fOg   = frac{2}({sWfg,sOfg,sGfg}, sTOfg, cfg    );
%     
%     fGc   = frac{3}({sWc ,sOc ,sGc }, sTc  , [c,c,c]);
%     fGv   = frac{3}({sWfv,sOfv,sGfv}, sTGfv, cfv    );
%     fGg   = frac{3}({sWfg,sOfg,sGfg}, sTGfg, cfg    );
    
    % Water equation + ncomp component equations
    [eqs, types, names] = deal(cell(1, 2*ncomp + model.water));
    mobOfg = mobO(sOfg,sTOfg,cfg(:,oind));
    mobGfg = mobG(sGfg,sTGfg,cfg(:,gind));

    % Water equation-------------------------------------------------------
    if model.water
        wIx = ncomp + 1;
        
        fWc   = frac{1}({sWc ,sOc ,sGc }, sTc  , [c,c,c]);
        
        % Cell values
        mobWc = mobW(sWc,sTc,c);
        mobOc = mobO(sOc,sTc,c);
        mobGc = mobO(sGc,sTc,c);
        % Face values
        fWfv   = frac{1}({sWfv,sOfv,sGfv}, sTWfv, cfv    );
        fWfg   = frac{1}({sWfg,sOfg,sGfg}, sTWfg, cfg    );
        mobWfg = mobW(sWfg,sTWfg,cfg(:,1));
        mobOfg = mobO(sOfg,sTOfg,cfg(:,2));
        mobGfg = mobG(sGfg,sTGfg,cfg(:,3));
        % Accumulation term
        acc = @(psi) (pvMult(c) .*rock.poro(c).*rhoW(c) .*sWc - ...
                      pvMult0(c).*rock.poro(c).*rhoW0(c).*sWc0).*psi/dt;
        % Convection term
        conv = @(gradPsi) rhoW(c).*fWc.*(disc.dot(vTc(c,:),gradPsi) ...
                        + mobOc.*disc.dot(TgWc(c,:) - TgOc(c,:),gradPsi) ...
                        + mobGc.*disc.dot(TgWc(c,:) - TgGc(c,:),gradPsi));
        integrand = @(psi, gradPsi) acc(psi) - conv(gradPsi);
        % Integrate integrand*psi{dofNo} over all cells for dofNo = 1:nDof
        cellIntegralW = disc.cellInt(integrand, [], state, sWdof);
        % Flux term
        integrand = @(psi) ...
            (sWfv.*sTWfv.*rhoW(cfv(:,1)).*fWfv.*vT(f) ...
                  + sWfg.*rhoW(cfg(:,1)).*fWfg.*mobOfg.*(TgW(f) - TgO(f)) ...
                  + sWfg.*rhoW(cfg(:,1)).*fWfg.*mobGfg.*(TgW(f) - TgG(f))).*psi;
        % Integrate integrand*psi{dofNo} over all cells surfaces for dofNo = 1:nDof
        faceIntegralW = disc.faceFluxInt(integrand, [], state, sWdof);
        % Sum integrals
        water = cellIntegralW + faceIntegralW;
        % Add well contributions
        if ~isempty(W)
            ix = disc.getDofIx(state, Inf, wc);
            water(ix) = water(ix) - srcWW(ix);
        end
        
        eqs{wIx} = water;
        names{wIx} = 'water';
        types{wIx} = 'cell';
        
    end
    
    %----------------------------------------------------------------------
    
    % Component equations--------------------------------------------------
    if model.water
        sc = {sWc ,sOc ,sGc };
        sfv = {sWfv,sOfv,sGfv};
        sfg = {sWfg,sOfg,sGfg};
    else
        sc = {sOc ,sGc };
        sfv = {sOfv,sGfv};
        sfg = {sOfg,sGfg};
    end
    fOc   = frac{oind}(sc, sTc  , [c,c,c]);
    fOfv   = frac{oind}(sfv, sTOfv, cfv    );
    fOfg   = frac{oind}(sfg, sTOfg, cfg    );
    
    fGc   = frac{gind}(sc, sTc  , [c,c,c]);
    fGfv   = frac{gind}(sfv, sTGfv, cfv    );
    fGfg   = frac{gind}(sfg, sTGfg, cfg    );
    
    mobOc = mobO(sOc,sTc,c);
    mobGc = mobO(sGc,sTc,c);
    
%     sOc = sOc./sTc;
%     xMc = cellfun(@(x) x.*sTc, 
    
    for cNo = 1:ncomp
        acc = @(psi) (pvMult(c) .*rock.poro(c).*rhoO(c) .*sOc .*xMc{cNo} - ...
                      pvMult0(c).*rock.poro(c).*rhoO0(c).*sOc0.*xMc0{cNo} + ...     
                      pvMult(c) .*rock.poro(c).*rhoG(c) .*sGc .*yMc{cNo}  - ...
                      pvMult0(c).*rock.poro(c).*rhoG0(c).*sGc0.*yMc0{cNo}).*psi/dt;
        if model.water
        conv = @(gradPsi) sTc.*rhoO(c).*xMc{cNo}.*fOc.*(...
                            disc.dot(vTc(c,:),gradPsi) ...
                          + mobWc.*disc.dot(TgOc(c,:) - TgWc(c,:),gradPsi) ...
                          + mobGc.*disc.dot(TgOc(c,:) - TgGc(c,:),gradPsi)) ...
                        + sTc.*rhoG(c).*yMc{cNo}.*fGc.*(...
                            disc.dot(vTc(c,:),gradPsi) ...
                          + mobWc.*disc.dot(TgGc(c,:) - TgWc(c,:),gradPsi) ...
                          + mobOc.*disc.dot(TgGc(c,:) - TgOc(c,:),gradPsi));
        else
        conv = @(gradPsi) sTc.*rhoO(c).*xMc{cNo}.*fOc.*(...
                            disc.dot(vTc(c,:),gradPsi) ...
                          + mobGc.*disc.dot(TgOc(c,:) - TgGc(c,:),gradPsi)) ...
                        + sTc.*rhoG(c).*yMc{cNo}.*fGc.*(...
                            disc.dot(vTc(c,:),gradPsi) ...
                          + mobOc.*disc.dot(TgGc(c,:) - TgOc(c,:),gradPsi));
        end
        
        integrand = @(psi, gradPsi) acc(psi) - conv(gradPsi);

        cellIntegral = disc.cellInt(integrand, [], state, sOdof);
        if model.water
        integrand = @(psi) ...
        ( ...
         sTOfv.*rhoO(cfv(:,2)).*xMOfv{cNo}.*fOfv.*vT(f) ...
              + rhoO(cfg(:,2)).*xMOfv{cNo}.*fOfg.*(mobWfg.*(TgO(f) - TgW(f)) ...
                                                 + mobGfg.*(TgO(f) - TgG(f))) + ...
         sTGfv.*rhoG(cfv(:,3)).*fGfv.*yMGfv{cNo}.*vT(f) ...
              + rhoG(cfg(:,3)).*fGfg.*yMGfv{cNo}.*(mobWfg.*(TgG(f) - TgW(f)) ...
                                                 + mobOfg.*(TgG(f) - TgO(f))) ...
        ).*psi;
        else
        integrand = @(psi) ...
        ( ...
         sTOfv.*(rhoO(cfv(:,oind)).*xMOfv{cNo}.*fOfv.*vT(f) ...
               + rhoO(cfg(:,oind)).*xMOfv{cNo}.*fOfg.*(mobGfg.*(TgO(f) - TgG(f)))) + ...
         sTGfv.*(rhoG(cfv(:,gind)).*fGfv.*yMGfv{cNo}.*vT(f) ...
               + rhoG(cfg(:,gind)).*fGfg.*yMGfv{cNo}.*(mobOfg.*(TgG(f) - TgO(f)))) ...
        ).*psi;
        end
        % Integrate integrand*psi{dofNo} over all cells surfaces for dofNo = 1:nDof
        faceIntegral = disc.faceFluxInt(integrand, [], state, sOdof);

        eqs{cNo} = cellIntegral + faceIntegral;
        if ~isempty(W)
            ix = disc.getDofIx(state, Inf, wc);
            eqs{cNo}(ix) = eqs{cNo}(ix) - srcW{cNo}(ix);
        end
        names{cNo} = compFluid.names{cNo};
        types{cNo} = 'cell';
    end
    
    %----------------------------------------------------------------------

if model.extraStateOutput
    bO = rhoO./fluid.rhoOS;
    bG = rhoG./fluid.rhoGS;
    state = model.storebfactors(state, bW, bO, bG);
%     state = model.storeMobilities(state, mobW, mobO, mobG);
end
state = model.storeDensities(state, rhoW, rhoO, rhoG);


% Fugacity
% tPhc = find(twoPhase);
zdof = expandMatrixToCell(zdof);
[zc] = deal(cell(size(xdof)));
% [~, xitPhc, tPhcNow] = disc.getCubature(tPhc, 'volume');
[zc{:}] = disc.evaluateDGVariable(xic, c, state, zdof{:});


% if ~isempty(twoPhaseIx) && ~opt.resOnly
%     offset = ncomp+model.water;
%     
%     [~, ~, ctPh] = disc.getCubature(find(twoPhase), 'volume');
% %     tp = find(twoPhase);
%     
%     
%     for i = 1:ncomp
%         cellJacMap{i + offset} = ctPh;%find(twoPhase);
%     end
%     
% % %     cellno = rldecode(a, b);
% %     cellNo = (1:G.cells.num);
% %     cellNo = cellNo(c);
%     cellNo = rldecode((1:G.cells.num)', state.nDof, 1);
%     
%     for i = 2:offset
%         cellJacMap{i} = cellNo;
%     end
% end
if 1
    [~, xic, c] = disc.getCubature((1:G.cells.num)', 'volume');
    [~,n] = rlencode(c);

    [f_L, f_V] = deal(cell(ncomp,1));
    [f_L{:}, f_V{:}] = deal(sO(c)*0);

    for xNo = 1:n(1)
        ind = xNo:n(1):size(xic,1);
        xx = xic(ind,:);
        cc = c(ind);

        [sOc , sGc , xc{:} , yc{:} , zc{:}] = disc.evaluateDGVariable(xx, cc, state , ...
                                 sOdof , sGdof , xdof{:} , ydof{:}, zdof{:});
%         xc = ensureMinimumFraction(xc);
%         yc = ensureMinimumFraction(yc);
        [xM,  yM,  rhoO,  rhoG,  muO,  muG, f_LL, f_VV, xM0, yM0, rhoO0, rhoG0] = ...
            model.getTimestepPropertiesEoS(state, state0, p, temp, xc, yc, zc, sOc, sGc, cellJacMap);
        for cNo = 1:ncomp
            f_L{cNo}(ind) = f_LL{cNo};
            f_V{cNo}(ind) = f_VV{cNo};
        end
    end
elseif 0
    % Compute properties and fugacity
    st = state;
    st.T = st.T(c);
    st.Z_L = st.Z_L(c);
    st.Z_V = st.Z_V(c);

    st0 = state0;
    st0.T = st0.T(c);
    st0.Z_L = st0.Z_L(c);
    st0.Z_V = st0.Z_V(c);
    st0.x = st0.x(c,:);
    st0.y = st0.y(c,:);
    st0.pressure = st0.pressure(c);

    [xM,  yM,  rhoO,  rhoG,  muO,  muG, f_L, f_V, xM0, yM0, rhoO0, rhoG0] = ...
        model.getTimestepPropertiesEoS(st, st0, p(c), temp(c), xc, yc, zc, sOc, sGc, cellJacMap);
end

z_tol = 1e-8;
for i = 1:ncomp
    cix = i + ncomp + model.water;
    names{cix}= ['f_', compFluid.names{i}];
    types{cix} = 'fugacity';
    
    integrand = @(psi, gradPsi) (f_L{i} - f_V{i}).*psi/barsa;
    
    fugacity = disc.cellInt(integrand, (1:G.cells.num)', state, sOdof);
    
    ix = disc.getDofIx(state, Inf, twoPhase);
    eqs{cix} = fugacity(ix)./rldecode(G.cells.volumes(twoPhase), state.nDof(twoPhase));
%     eqs{cix} = (f_L{i}(twoPhase) - f_V{i}(twoPhase))/barsa;
    absent = state.components(twoPhase, i) <= 10*z_tol;
    if model.water
        absent = absent | pureWater(twoPhase);
    end
    if any(absent) && isa(eqs{cix}, 'ADI')
        absIx = disc.getDofIx(state, Inf, absent);
        eqs{cix}.val(absIx) = 0;
    end    
end
massT = model.getComponentScaling(state0);
scale = (dt./op.pv)./massT;
if model.water
    wscale = dt./(op.pv*mean(double(rhoW0)));
    wscale = rldecode(wscale, state.nDof);
    eqs{wIx} = eqs{wIx}.*wscale;
end
scale = rldecode(scale, state.nDof);
for cNo = 1:ncomp
    eqs{cNo} = eqs{cNo}.*scale;
end


if model.reduceLinearSystem
    problem = ReducedLinearizedSystem(eqs, types, names, primaryVars, state, dt);
    
    problem.keepNum = sum(state.nDof).*(ncomp + model.water);
%     problem.keepNum = model.G.cells.num*(ncomp+model.water);
    
    problem.reorder = reorder;
else
    problem = LinearizedProblem(eqs, types, names, primaryVars, state, dt);
end

if 0
    
    eqsDG = eqs;
    load('fv.mat');
    [d, djac] = deal(cell(size(eqs)));
    
    for eqNo = 1:numel(eqs)
        d{eqNo} = eqs{eqNo} - eqsDG{eqNo};
        djac{eqNo} = cellfun(@(d) norm(d, 'fro'), d{eqNo}.jac);
        djacS{eqNo} = cellfun(@(d,e) norm(d, 'fro')./norm(e, 'fro'), d{eqNo}.jac, eqs{eqNo}.jac);
    end
    dval = cellfun(@(d) norm(d.val), d);
    
    
%     djac = cellfun(@(d) norm(d.val), d, ');
    
end

if 0
    ix = disc.getDofIx(state, 2, Inf);
    for eqNo = 1:4
        v(:,eqNo) = eqs{eqNo}(ix).val;
    end
    plotToolbar(G, v, 'plot1d', true);
end
% [sum(pureLiquid), sum(pureVapor), sum(twoPhase)]
end

% Expand single scalar values to one per cell------------------------------
function v = expandSingleValue(v,G)
    if numel(double(v)) == 1
        v = v*ones(G.cells.num,1);
    end
end
%--------------------------------------------------------------------------

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
