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

    % Prepare state for simulation-----------------------------------------
    if opt.iteration == 1 && ~opt.resOnly 
        if model.tryMaxDegree
            % If we are at the first iteration, we try to solve using
            % maximum degree in all cells
            state.degree(~G.cells.ghost) = disc.degree;
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

    fluid = model.fluid;
    compFluid = model.EOSModel.fluid;

    % Properties at current timestep
    [p, sWdof, sOdof, sGdof, xdof, ydof, temp, wellSol] = model.getProps(state, ...
        'pressure', 'swdof', 'sodof', 'sgdof', 'xdof', 'ydof', 'T', 'wellSol');
    assert(all(p>0), 'Pressure must be positive for compositional model');

    [p0, sWdof0, sOdof0, sGdof0, xdof0, ydof0, temp0, wellSol0] = model.getProps(state0, ...
        'pressure', 'swdof', 'sodof', 'sgdof', 'xdof', 'ydof', 'T', 'wellSol');

    [pureLiquid, pureVapor, twoPhase] = model.getFlag(state);
    pureLiquidIx = disc.getDofIx(state, Inf, pureLiquid);
    pureVaporIx  = disc.getDofIx(state, Inf, pureVapor);
    
    if 1
        % TODO: Do the same for dofs
        sO = state.s(:,2);
        sG = state.s(:,3);
        stol = 1e-8;
        pureWater = sO + sG < stol;
        sO(~pureVapor & pureWater) = stol;
        sG(~pureLiquid & pureWater) = stol;

        sO0 = state0.s(:,2);
        sG0 = state0.s(:,3);
        [pureLiquid0, pureVapor0, twoPhase0] = model.getFlag(state0);
        pureWater0 = sO0 + sG0 < stol;
        sO0(~pureVapor0 & pureWater0) = stol;
        sG0(~pureLiquid0 & pureWater0) = stol;
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
    if 0
        % TODO: Make function ensureMinimumFraction for dofs
        x = ensureMinimumFraction(x);
        [y, z_tol] = ensureMinimumFraction(y);
    end
    xdof = expandMatrixToCell(xdof);
    ydof = expandMatrixToCell(ydof);

    ncomp = model.EOSModel.fluid.getNumberOfComponents();
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
    sodof = sOdof(twoPhase);

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
        [X, xdof{1:ncomp-1}, wdof{1:  ncomp-1}, sodof] = model.AutoDiffBackend.initVariablesAD(...
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
ix = disc.getDofIx(state, 1);

[xdof{end}, wdof{end}] = deal(zeros(size(double(xdof{1}))));
xdof{end}(ix) = 1;
wdof{end}(ix) = 1;

% xdof{end} = 1;
% wdof{end} = 1;
for i = 1:ncomp-1
    xdof{end} = xdof{end} - xdof{i};
    if any(twoPhase)
        wdof{end} = wdof{end}-wdof{i};
    end
end

for i = 1:ncomp
    ydof{i} = xdof{i};
    if any(twoPhase)
        ydof{i}(twoPhaseIx) = wdof{i};
    end
    xdof{i}(pureVaporIx) = double(xdof{i}(pureVaporIx));
    ydof{i}(pureLiquidIx) = double(ydof{i}(pureLiquidIx));
end

cellJacMap = cell(numel(primaryVars), 1);
% toReorder = 1:nc;

if isempty(twoPhaseIx) || opt.resOnly
    reorder = [];
else
    % TODO: Implement for dg. Remember twoPhaseIx not the same as in
    % original code
    n2ph = nnz(twoPhaseIx);
    nVars = sum(sample.getNumVars());
    reorder = 1:nVars;
    start = twoPhaseIx + nc;
    stop = (nVars-n2ph+1):nVars;
    
    reorder(start) = stop;
    reorder(stop) = start;
    
    offset = ncomp+model.water;
    for i = 1:ncomp
        cellJacMap{i + offset} = twoPhaseIx;
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
if model.water
    sW = model.disc.getCellMean(sWdof, state);
end
sO = model.disc.getCellMean(sOdof, state);
sG = model.disc.getCellMean(sGdof, state);
sT = model.disc.getCellMean(sTdof, state);

% Compute properties and fugacity
[xM,  yM,  rhoO,  rhoG,  muO,  muG, f_L, f_V, xM0, yM0, rhoO0, rhoG0] = ...
                  model.getTimestepPropertiesEoS(state, state0, p, temp, xdof, ydof, zdof, sOdof, sGdof, cellJacMap);

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
    mobW  = @(sW,sT,c) mobMult(c).*fluid.krW(sW./sT)./muG(c);

    gW    = rhoWf.*gdz;

    if isfield(fluid, 'pcOW')
        gW = gW + op.Grad(fluid.pcOW(sW));
    end
    gW = P*gW;
    sWT       = sW.*sT;
    g         = {gW, gO, gG};
    mob       = {mobW, mobO, mobG};
    rho       = {rhoW, rhoO, rhoG};
    pressures = {pW, p, p};
    sdof      = {sWdof, sOdof, sGdof};
else
    [rhoW, rhoW0, mobW, bW, sWT] = deal([]);
    g         = {gO, gG};
    mob       = {mobO, mobG};
    rho       = {rhoO, rhoG};
    sdof      = {sOdof, sGdof};
    pressures = {p, p};
end

    % Evaluate saturation/compositions at cubature points------------------
    % Cell cubature points
    xdof0 = expandMatrixToCell(xdof0);
    ydof0 = expandMatrixToCell(ydof0);
    [~, xic, c] = disc.getCubature((1:G.cells.num)', 'volume');
    [xc, yc, xc0, yc0, xOfv, xGfv, yOfv, yGfv] = deal(cell(size(xdof)));
    [sWc , sOc , sGc , sTc , xc{:} , yc{:} ] = disc.evaluateDGVariable(xic, c, state , ...
                             sWdof , sOdof , sGdof , sTdof , xdof{:} , ydof{:});
    [sWc0, sOc0, sGc0, sTc0, xc0{:}, yc0{:}] = disc.evaluateDGVariable(xic, c, state0, ...
                             sWdof0, sOdof0, sGdof0, sTdof0, xdof0{:}, ydof0{:});
    
    % Face cubature points
    [~, xif, ~, f] = disc.getCubature((1:G.cells.num)', 'surface');
    % Upstream cells
    [~, ~, cfv, cfg] = disc.getSaturationUpwind(f, xif, T, flux, g, mob, sdof, state);
    % Water saturation
    [sWfv , sTWfv] = disc.evaluateDGVariable(xif, cfv(:,1), state, sWdof, sTdof, xdof{:}, ydof{:});
    [sWfg , sTWfg] = disc.evaluateDGVariable(xif, cfg(:,1), state, sWdof, sTdof);
    
    
    % Oil saturation
    [sOfV , sTOfV, xOfv{:}, yOfv{:}] = disc.evaluateDGVariable(xif, cfv(:,2), state, sOdof, sTdof, xdof{:}, ydof{:});
    [sOfG , sTOfg] = disc.evaluateDGVariable(xif, cfg(:,2), state, sOdof, sTdof);
    % Gas saturation
    [sGfV , sTGfV, xGfv{:}, yGfv{:}] = disc.evaluateDGVariable(xif, cfv(:,3), state, sGdof, sTdof, xdof{:}, ydof{:});
    [sGfG , sTGfg] = disc.evaluateDGVariable(xif, cfg(:,3), state, sGdof, sTdof);
    
    
    

% [xM, yM] = deal(cell(ncom),1);
% for cNo = 1:ncomp
%     
% end

% for i = 1:ncomp
%     xM{i} = xM{i}.*sT;
%     yM{i} = yM{i}.*sT;
%     
%     xM0{i} = xM0{i}.*sT0;
%     yM0{i} = yM0{i}.*sT0;
% end

if model.extraStateOutput
    bO = rhoO./fluid.rhoOS;
    bG = rhoG./fluid.rhoGS;
    state = model.storebfactors(state, bW, bO, bG);
%     state = model.storeMobilities(state, mobW, mobO, mobG);
end
state = model.storeDensities(state, rhoW, rhoO, rhoG);


components = getComponentsTwoPhaseSimpleWater(model, rho, sWT, xM, yM);

upstr = model.operators.faceUpstr;
[q_phase, q_components] = computeSequentialFluxes(...
    state, g, vT, T, mob, rho, components, upstr, model.upwindType);

pv = pvMult.*model.operators.pv;
pv0 = pvMult0.*model.operators.pv;

compFlux = zeros(size(model.operators.N, 1), ncomp);

% water equation + n component equations
[eqs, types, names] = deal(cell(1, 2*ncomp + model.water));
for i = 1:ncomp
    names{i} = compFluid.names{i};
    types{i} = 'cell';
    vi = q_components{i};
    eqs{i} = (1/dt).*( ...
                    pv.*rhoO.*sO.*xM{i} - pv0.*rhoO0.*sO0.*xM0{i} + ...
                    pv.*rhoG.*sG.*yM{i} - pv0.*rhoG0.*sG0.*yM0{i}) + s.Div(vi);
      
   compFlux(:, i) = double(vi);
end
state.componentFluxes = compFlux;

if model.water
    wix = ncomp+1;
    vw = q_components{ncomp+1};
    eqs{wix} = (1/dt).*(pv.*rhoW.*sW.*sT - pv0.*rhoW0.*sW0.*sT0) + s.Div(vw);
    names{wix} = 'water';
    types{wix} = 'cell';
    state = model.storeFluxes(state, q_phase{:});
else
    state = model.storeFluxes(state, [], q_phase{:});
end

comps = cellfun(@(x, y) {x, y}, xM, yM, 'UniformOutput', false);

[eqs, state] = model.addBoundaryConditionsAndSources(eqs, names, types, state, ...
                                                 pressures, sat, mob, rho, ...
                                                 {}, comps, ...
                                                 drivingForces);


% Finally, add in and setup well equations
if ~isempty(W)
    mflux = sum(vertcat(wellSol.components), 2);
    wflux = sum(vertcat(wellSol.flux), 2);

    perf2well = getPerforationToWellMapping(W);
    wc    = vertcat(W.cells);
    w_comp = vertcat(W.components);
    a = w_comp(perf2well, :).*repmat(compFluid.molarMass, numel(wc), 1);
    w_comp = bsxfun(@rdivide, a, sum(a, 2));

    isInj = wflux > 0;
    compWell = vertcat(W.compi);
    compPerf = compWell(perf2well, :);

    x_comp = cellfun(@(v) v(wc), xM, 'UniformOutput', false);
    y_comp = cellfun(@(v) v(wc), yM, 'UniformOutput', false);

    mobOw = mobO(wc);
    mobGw = mobG(wc);

    rhoOw = rhoO(wc);
    rhoGw = rhoG(wc);

    if model.water
        mobWw = mobW(wc);
        rhoWw = rhoW(wc);
        totMobw = mobWw + mobOw + mobGw;
        f_w_w = mobWw./totMobw;
        f_w_w = sT(wc).*f_w_w;

        f_w_w(isInj) = compPerf(isInj, 1);
        rWqW = rhoWw.*f_w_w.*wflux;
        eqs{wix}(wc) = eqs{wix}(wc) - rWqW;
    else
        totMobw = mobOw + mobGw;
    end
    f_o_w = mobOw./totMobw;
    f_g_w = mobGw./totMobw;
    rOqO = rhoOw.*f_o_w.*wflux;
    rGqG = rhoGw.*f_g_w.*wflux;


    rOqO(isInj) = compPerf(isInj, 1 + model.water).*mflux(isInj);
    rGqG(isInj) = compPerf(isInj, 2 + model.water).*mflux(isInj);

    sources = cell(ncomp, 1);
    compSrc = zeros(numel(wc), ncomp);
    for i = 1:ncomp
        src =       (rOqO + rGqG).*w_comp(:, i).*isInj...
                   +(rOqO.*x_comp{i} + rGqG.*y_comp{i}).*~isInj;
        eqs{i}(wc) = eqs{i}(wc) - src;
        compSrc(:, i) = double(src);
        sources{i} = src;    
    end

    qg = double(rGqG)./fluid.rhoGS;
    qo = double(rOqO)./fluid.rhoOS;

    isProd = ~isInj;
    qo(isProd) = qo(isProd).*sT(wc(isProd));
    qg(isProd) = qg(isProd).*sT(wc(isProd));
    if model.water
        qw = double(rWqW)./fluid.rhoWS;
        qw(isProd) = qw(isProd).*sT(wc(isProd));
    end
    for i = 1:numel(W)
        state.wellSol(i).components = compSrc(perf2well == i, :);
        state.wellSol(i).qGs = sum(qg(perf2well == i));
        state.wellSol(i).qOs = sum(qo(perf2well == i));
        if model.water
            state.wellSol(i).qWs = sum(qw(perf2well == i));
        end

        state.wellSol(i).qTr = sum(wflux(perf2well == i));
        state.wellSol(i).status = true;
        state.wellSol(i).type = W(i).type;
        state.wellSol(i).val = W(i).val;
    end
end

for i = 1:ncomp
    ix = i + ncomp + model.water;
    names{ix}= ['f_', compFluid.names{i}];
    types{ix} = 'fugacity';
    eqs{ix} = (f_L{i}(twoPhase) - f_V{i}(twoPhase))/barsa;
    absent = state.components(twoPhase, i) <= 10*z_tol;
    if model.water
        absent = absent | pureWater(twoPhase);
    end
    if any(absent) && isa(eqs{ix}, 'ADI')
        eqs{ix}.val(absent) = 0;
    end    
end
massT = model.getComponentScaling(state0);
scale = (dt./s.pv)./massT;
if model.water
    wscale = dt./(s.pv*mean(double(rhoW0)));
    eqs{wix} = eqs{wix}.*wscale;
end

for i = 1:ncomp
    eqs{i} = eqs{i}.*scale;
end


if model.reduceLinearSystem
    problem = ReducedLinearizedSystem(eqs, types, names, primaryVars, state, dt);
    problem.keepNum = model.G.cells.num*(ncomp+model.water);
    problem.reorder = reorder;
else
    problem = LinearizedProblem(eqs, types, names, primaryVars, state, dt);
end
end

% Expang single scalar values to one per cell------------------------------
function v = expandSingleValue(v,G)
    if numel(double(v)) == 1
        v = v*ones(G.cells.num,1);
    end
end
%--------------------------------------------------------------------------

%{
Copyright 2009-2018 SINTEF ICT, Applied Mathematics.

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
