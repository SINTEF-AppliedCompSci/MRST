function [dx, dy, ds, dL, twoPhase, w] = getUpdatesFromPressure(model, state, dp)
    [~, ~, twoPhase] = model.getFlag(state);
    noHydrocarbons = state.s(twoPhase, 1+model.water) == 0 & state.s(twoPhase, 2+model.water) == 0;
    
    if ~any(twoPhase)
%         [dx, dy, ds, dL] = deal([]);
        dL = [];
        w = [];
        ncomp = size(state.components, 2);
        nph = size(state.s, 2);
        dx = zeros(0, ncomp);
        dy = zeros(0, ncomp);
        ds = zeros(0, nph);
        return
    end
    
    x0 = state.x(twoPhase, :);
    y0 = state.y(twoPhase, :);
%     z = expandMatrixToCell(state.components, twoPhase);
    z0 = state.z0(twoPhase, :);

    x = expandMatrixToCell(x0);
    y = expandMatrixToCell(y0);
%     z = expandMatrixToCell(state.components, twoPhase);
    z = expandMatrixToCell(z0);

%     state.s = state.s./sum(state.s, 2);
    ncomp = numel(x);
    
    dp = dp(twoPhase);
    temp = state.T(twoPhase);
    p = state.pressure(twoPhase);

    
    Z_L = state.Z_L(twoPhase);
    Z_V = state.Z_V(twoPhase);
    L0 = state.L(twoPhase);
    sO = state.s(twoPhase, 1+model.water);
    sG = state.s(twoPhase, 2+model.water);
    if model.water
        sW = state.s(twoPhase, 1);
    else
        sW = 0;
    end
    
%     [p, x{1:end}, y{1:end}, sO, sG, L] = initVariablesADI(p, x{1:end}, y{1:end}, sO, sG, L0);
    includeWater = model.water;
    if includeWater
        [p, x{1:end}, y{1:end}, sW, sO, sG] = initVariablesADI(p, x{1:end}, y{1:end}, sW, sO, sG);
    else
        [p, x{1:end}, y{1:end}, sO, sG] = initVariablesADI(p, x{1:end}, y{1:end}, sO, sG);
    end
    
    eos = model.EOSModel;
    
    tmp = struct('Z_L', Z_L, 'Z_V', Z_V);
    [Z_L, Z_V, f_L, f_V] = eos.getProperties(p, temp, x, y, z, sO, sG, tmp);
    rhoO_m = model.PropertyModel.computeMolarDensity(p, x, Z_L, temp, true);
    rhoG_m = model.PropertyModel.computeMolarDensity(p, y, Z_V, temp, false);
    
    rhoO = model.PropertyModel.computeDensity(p, x, Z_L, temp, true);
    rhoG = model.PropertyModel.computeDensity(p, y, Z_V, temp, false);

    mol_L = sO.*rhoO_m;
    mol_G = sG.*rhoG_m;
    mol_T = mol_G + mol_L;
    
    eqs = cell(2*ncomp + 3, 1);
    eqs{ncomp*2 + 1} = 0;
    eqs{2*ncomp} = 1;
    eqs{2*ncomp+1} = 1;
    for i = 1:ncomp
        eqs{i} = (f_L{i} - f_V{i})/barsa;
        if i < ncomp
            eqs{i+ncomp} = mol_L.*x{i} + mol_G.*y{i} - mol_T.*z{i};
        end
        eqs{2*ncomp} = eqs{2*ncomp} - x{i};
        eqs{2*ncomp+1} = eqs{2*ncomp+1} - y{i};
    end
    
    if includeWater
        rhoW = model.fluid.bW(p).*model.fluid.rhoWS;
        eqs{2*ncomp + 2} = sW.*rhoW./(sW.*rhoW + sO.*rhoO + sG.*rhoG) - state.sM(twoPhase);
    end
    eqs{2*ncomp + 2 + includeWater} = sO + sG + sW - 1;
    
    eqn = combineEquations(eqs{:});

    n = nnz(twoPhase);
    
    A = -eqn.jac{1};
    Ap = A(:, 1:n);
    Ar = A(:, n+1:end);
    b = eqn.val;
    
    r = (b - Ap*dp);
    res = Ar\r;
    
    phaseix = 1:n*ncomp;
    
    dx = res(phaseix);
    dy = res(phaseix + n*ncomp);
    
    dx = reshape(dx, [], ncomp);
    dy = reshape(dy, [], ncomp);
    
    nPh = (2 + includeWater);
    nsat = nPh*n;
    ds = res((2*n*ncomp+1):(2*n*ncomp+nsat));
    ds = reshape(ds, [], nPh);
    
    dL = [];
    
    unit = @(x) min(max(x, 0), 1);
    x = unit(x0 + dx);
    y = unit(y0 + dy);
    s0 = state.s(twoPhase, :);
    s = max(min(s0 + ds, -0.01), 1.01);
    
    wx = (x - x0)./dx;
    wy = (y - y0)./dy;
    ws = (s - s0)./ds;
    
    wx(abs(dx) < 0.1*model.nonlinearTolerance) = 1;
    wy(abs(dy) < 0.1*model.nonlinearTolerance) = 1;
    ws(abs(ds) < 0.1*model.nonlinearTolerance) = 1;
    ws = min(ws, [], 2);
    
    w = min(wx, wy);
    w = min(w, [], 2);
%     w = min(w, ws);
    
    dx = bsxfun(@times, dx, w);
    dy = bsxfun(@times, dy, w);
    ds = bsxfun(@times, ds, w);
    
    dx(noHydrocarbons, :) = 0;
    dy(noHydrocarbons, :) = 0;
    
%     dL = res((2*n*ncomp+nsat+1):end);
    
%     L = L0 + dL;
%     V = 1-L;
%     L = (z0 - y0 - dy)./(x0 + dx - y0 - dy);
%     
%     [(x0 + dx).*(L0 + dL) + (y0 + dy).*(1 - L0 - dL)]./z0
%     
%     [(x0 + w.*dx).*(L0 + w.*dL) + (y0 + w.*dy).*(1 - L0 - w.*dL)]./z0
end
