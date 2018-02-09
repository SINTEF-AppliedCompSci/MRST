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
    
%     acc = 1/dt*op.M.*(sWdof - sWdof0);
%     
%     flux = 0;

    [k, nDof] = dgBasis(model.degree, model.G.griddim);

    sW  = cell(model.G.cells.num,1);
%     sW0 = @(cell,x) 0*x;
    
    w = 1;
    psi = @(x) [w*prod(x.^[0,0],2), w*prod(x.^[1,0],2), w*prod(x.^[0,1],2)];
    
    psi = {@(x) w*prod(x.^[0,0],2), @(x) w*prod(x.^[1,0],2), @(x) w*prod(x.^[0,1],2)};
    
    sW = satfun(sWdof, psi, model.G);
    ss = sW([1,1]);

    for cNo = 1:model.G.cells.num
        sWc = 0;
        for dofNo = 1:nDof
            psi = Polynomial(k(dofNo,:));
            sWc  = sWc + sWdof((cNo-1)*dofNo + dofNo)*psi;
        end
        sW{cNo} = sWc;
        sloc = zeros(nDof,1);
        [intFun, w, x] = makeCellIntegrator(model.G, cNo, model.degree, 'tri');
        for dofNo = 1:nDof
            psi = Polynomial(k(dofNo,:));
            sloc(dofNo) = intFun(sWc*psi);
        end
    end
    
    




%     sWdof.*
%     psi(x);
    
    sO = 1-sW;
    
    [krW, krO] = model.evaluateRelPerm({sW, sO});
    
    [pvMult, transMult, mobMult, pvMult0] = getMultipliers(model.fluid, p, p0);

    % Modifiy relperm by mobility multiplier (if any)
    krW = mobMult.*krW; krO = mobMult.*krO;

    % Compute transmissibility
    T = s.T.*transMult;
    
    [vW, bW, mobW, rhoW, pW, upcw, dpW] = getFluxAndPropsWater_BO(model, p, sW, krW, T, gdz);
    
    for dofNo = 1:nDof
        psi = Polynomial(k(dofNo,:));
        ind = (1:nDof:G.cells.num) + dofNo - 1;
        accW(ind)  = op.cellInt(sW.*bW.*psi);
        fluxW(ind) = op.cellInt(fW(sW).*op.div(v*psi)) + ...
                     op.faceCellInt(fWhat(bWvW)*psi);
        
    end
    
    water = accW + fluxW;
    
%     accW = 1/dt.*op.M*(sWdof - sWdof0);


%     fluxW = op.cellInt(
%     
    % -------------------------------------------------------------------------
    sO = 1 - sW;
    [krW, krO] = model.evaluateRelPerm({sW, sO});

    % Multipliers for properties
    [pvMult, transMult, mobMult, pvMult0] = getMultipliers(model.fluid, p, p0);

    % Modifiy relperm by mobility multiplier (if any)
    krW = mobMult.*krW; krO = mobMult.*krO;

    % Compute transmissibility
    T = s.T.*transMult;

    % Gravity gradient per face
    gdz = model.getGravityGradient();

    % Evaluate water properties
    [vW, bW, mobW, rhoW, pW, upcw, dpW] = getFluxAndPropsWater_BO(model, p, sW, krW, T, gdz);

    % Evaluate oil properties
    [vO, bO, mobO, rhoO, pO, upco, dpO] = getFluxAndPropsOil_BO(model, p, sO, krO, T, gdz);

    gp = s.Grad(p);
    Gw = gp - dpW;
    Go = gp - dpO;

    if model.extraStateOutput
        state = model.storebfactors(state, bW, bO, []);
        state = model.storeMobilities(state, mobW, mobO, []);
    end

if ~isempty(W)
    wflux = sum(vertcat(wellSol.flux), 2);
    perf2well = getPerforationToWellMapping(W);
    wc = vertcat(W.cells);
    
    mobWw = mobW(wc);
    mobOw = mobO(wc);
    totMobw = mobWw + mobOw;

    f_w_w = mobWw./totMobw;
    f_o_w = mobOw./totMobw;

    isInj = wflux > 0;
    compWell = vertcat(W.compi);
    compPerf = compWell(perf2well, :);

    f_w_w(isInj) = compPerf(isInj, 1);
    f_o_w(isInj) = compPerf(isInj, 2);

    bWqW = bW(wc).*f_w_w.*wflux;
    bOqO = bO(wc).*f_o_w.*wflux;

    % Store well fluxes
    wflux_O = double(bOqO);
    wflux_W = double(bWqW);
    
    for i = 1:numel(W)
        perfind = perf2well == i;
        state.wellSol(i).qOs = sum(wflux_O(perfind));
        state.wellSol(i).qWs = sum(wflux_W(perfind));
    end

end

% Get total flux from state
flux = sum(state.flux, 2);
vT = flux(model.operators.internalConn);

% Stored upstream indices
[flag_v, flag_g] = getSaturationUpwind(model.upwindType, state, {Gw, Go}, vT, s.T, {mobW, mobO}, s.faceUpstr);
flag = flag_v;

upcw  = flag(:, 1);
upco  = flag(:, 2);

upcw_g = flag_g(:, 1);
upco_g = flag_g(:, 2);

mobOf = s.faceUpstr(upco, mobO);
mobWf = s.faceUpstr(upcw, mobW);

totMob = (mobOf + mobWf);
    
mobWf_G = s.faceUpstr(upcw_g, mobW);
mobOf_G = s.faceUpstr(upco_g, mobO);
mobTf_G = mobWf_G + mobOf_G;
f_g = mobWf_G.*mobOf_G./mobTf_G;
if opt.solveForWater
    f_w = mobWf./totMob;
    bWvW   = s.faceUpstr(upcw, bW).*f_w.*vT + s.faceUpstr(upcw_g, bO).*f_g.*s.T.*(Gw - Go);

    wat = (s.pv/dt).*(pvMult.*bW.*sW - pvMult0.*f.bW(p0).*sW0) + s.Div(bWvW);
    if ~isempty(W)
        wat(wc) = wat(wc) - bWqW;
    end

    eqs = {wat};
    names = {'water'};    
else
    f_o = mobOf./totMob;
    bOvO   = s.faceUpstr(upco, bO).*f_o.*vT + s.faceUpstr(upco_g, bO).*f_g.*s.T.*(Go - Gw);

    oil = (s.pv/dt).*( pvMult.*bO.*(1-sW) - pvMult0.*f.bO(p0).*(1-sW0) ) + s.Div(bOvO);
    if ~isempty(W)
        oil(wc) = oil(wc) - bOqO;
    end
    eqs = {oil};
    names = {'oil'};
end
types = {'cell'};
rho = {rhoW, rhoO};
mob = {mobW, mobO};
sat = {sW, sO};
[eqs, ~, src] = addBoundaryConditionsAndSources(model, eqs, names, types, state, ...
                                     {pFlow, pFlow}, sat, mob, rho, ...
                                     {}, {}, ...
                                     drivingForces);

if ~model.useCNVConvergence
    eqs{1} = eqs{1}.*(dt./s.pv);
end

problem = LinearizedProblem(eqs, types, names, primaryVars, state, dt);
end

function fun = satfun(sdof, psi, G)
    
    fun = @(x) x(:,1)*0;
    nDof = size(psi,2);
    for dofNo = 1:size(psi,2)
        ix = (1:nDof:G.cells.num*nDof) + dofNo - 1;
        fun = @(x) fun(x) + sdof(ix).*psi{dofNo}(x);
    end
    
end