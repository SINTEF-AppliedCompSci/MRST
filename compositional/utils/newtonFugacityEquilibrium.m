function [x, y, K, Z_L, Z_V, L, values] = newtonFugacityEquilibrium(model, P, T, z, K, L)
% Single step of the newton update for flash equations
%
% SYNOPSIS:
%   [x, y, K, Z_L, Z_V, L, values] = newtonFugacityEquilibrium(model, P, T, z, K, L)
%
% DESCRIPTION:
%   Perform a single Newton update 
%
% PARAMETERS:
%   eos - EquationOfState derived class instance.
%   p   - Pressures as a column vector. One value per cell.
%   T   - Temperatures as a column vector. One value per cell.
%   z   - Composition as a matrix with number of rows equal to the number
%         of components.
%   K   - Equilibrium constants as a matrix with number of rows equal to the number
%         of components.
%   L   - Column vector of liquid mole fractions. One value per cell.
%
% SEE ALSO:
%   `EquationOfStateModel`

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
    ncomp = model.getNumberOfComponents();
    ncell = numel(L);

    z = ensureMinimumFraction(z);
    
    x = model.computeLiquid(L, K, z);
    y = model.computeVapor(L, K, z);
    x = expandMatrixToCell(x);
    y = expandMatrixToCell(y);
    z = expandMatrixToCell(z);
    
    [xAD, yAD] = deal(cell(1, ncomp));
    [xAD{:}, yAD{:}, L_ad] = model.AutoDiffBackend.initVariablesAD(x{:}, y{:}, L);
    [eqs, f_L, f_V, Z_L, Z_V] = model.equationsEquilibrium(P, T, xAD, yAD, z, L_ad, [], []);

    eq = combineEquations(eqs{:});
    
    J = eq.jac{1};
    r = eq.val;
    upd = -J\r;
    cap = @(x) min(max(x, 0), 1);

    [dx, dy] = deal(cell(1, ncomp));
    
    dL = upd(end-ncell + 1: end);
    bad = isinf(dL);
    dL(bad) = sign(dL(bad));
    
    w = getRelax(L, dL, cap);
    for i = 1:ncomp
        subs = (ncell*(i-1) + 1):(ncell*i);
        dx{i} = upd(subs);
        dy{i} = upd(subs + ncomp*ncell);
        
        dx{i}(L == 0) = 0;
        dy{i}(L == 1) = 0;
        
        w_x = getRelax(x{i}, dx{i}, cap);
        w_y = getRelax(y{i}, dy{i}, cap);
        w = min(w, min(w_x, w_y));        
    end
    
    f_r = cell(1, ncomp);
    L = L + w.*dL;

    pure = L == 1 | L == 0;

    K = cell(1, ncomp);
    for i = 1:ncomp
        zmissing = z{i} == 0;
        x{i} = x{i} + w.*dx{i};
        y{i} = y{i} + w.*dy{i};
        x{i}(pure) = z{i}(pure);
        y{i}(pure) = z{i}(pure);
        
        f_r{i} = value(f_L{i})./value(f_V{i});
        f_r{i}(pure | zmissing) = 1;
        
        K{i} = y{i}./x{i};
        
        xmissing = x{i} == 0;
        gone = zmissing | xmissing | pure;
        K{i}(gone) = 1;
    end
    fn = @(v) cellfun(@(x) all(isfinite(x)), v);
    assert(all(fn(K)));
    assert(all(fn(x)));
    assert(all(fn(y)));
    assert(all(isfinite(L)));
    
    K = [K{:}];
    x = [x{:}];
    y = [y{:}];
    
    values = abs([f_r{:}] - 1);
end


function w = getRelax(x, dx, cap)
    w = (cap(x + dx) - x)./(dx);
    w(isnan(w)) = 0;
end

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

