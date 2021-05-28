function [dx, dy, ds, dL, twoPhase, w] = getUpdatesFromPressure(model, state, dp)
%Undocumented Utility Function

%{
Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.

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

    Z_L = value(state.Z_L);
    Z_V = value(state.Z_V);
    Z_L = Z_L(twoPhase);
    Z_V = Z_V(twoPhase);
    sO = state.s(twoPhase, 1+model.water);
    sG = state.s(twoPhase, 2+model.water);
    if model.water
        sW = state.s(twoPhase, 1);
    else
        sW = 0;
    end
    includeWater = model.water;
    if includeWater
        [p, x{1:end}, y{1:end}, sW, sO, sG] = model.AutoDiffBackend.initVariablesAD(p, x{1:end}, y{1:end}, sW, sO, sG);
    else
        [p, x{1:end}, y{1:end}, sO, sG] = model.AutoDiffBackend.initVariablesAD(p, x{1:end}, y{1:end}, sO, sG);
    end
    eos = model.EOSModel;
    [Z_L, Z_V, f_L, f_V] = eos.getCompressibilityAndFugacity(p, temp, x, y, z, Z_L, Z_V);
    rhoO_m = model.PropertyModel.computeMolarDensity(eos, p, x, Z_L, temp, true);
    rhoG_m = model.PropertyModel.computeMolarDensity(eos, p, y, Z_V, temp, false);
    
    rhoO = model.PropertyModel.computeDensity(eos, p, x, Z_L, temp, true);
    rhoG = model.PropertyModel.computeDensity(eos, p, y, Z_V, temp, false);

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
%     s0 = state.s(twoPhase, :);
%     
%     s_tol = 1e-3;
%     s = min(max(s0 + ds, -s_tol), 1 + s_tol);
    
    wx = abs(x - x0)./abs(dx);
    wy = abs(y - y0)./abs(dy);
%     ws = abs(s - s0)./abs(ds);
    
    wx(abs(dx) < 0.1*model.nonlinearTolerance) = 1;
    wy(abs(dy) < 0.1*model.nonlinearTolerance) = 1;
%     ws(abs(ds) < 0.1*model.nonlinearTolerance) = 1;
%     ws = min(ws, [], 2);
    
    w = min(wx, wy);
    w = min(w, [], 2);
%     w = min(w, ws);

    dx = bsxfun(@times, dx, w);
    dy = bsxfun(@times, dy, w);
    ds = bsxfun(@times, ds, w);
    
    dx(noHydrocarbons, :) = 0;
    dy(noHydrocarbons, :) = 0;
end
