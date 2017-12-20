function [stable, x, y] = phaseStabilityTest(eos, z, p, T, x, y, varargin)
% Perform a phase stability test for a mixture
%
% SYNOPSIS:
%   [stable, x, y] = phaseStabilityTest(eos, z, p, T, x, y)
%
% DESCRIPTION:
%   Perform a phase stability test for a multicomponent 
%
% PARAMETERS:
%   eos - EquationOfState derived class instance.
%   z   - Composition as a matrix with number of rows equal to the number
%         of components.
%   p   - Pressures as a column vector
%   T   - Temperatures as a column vector
%   x   - Initial guess for liquid mole fractions as a matrix with number
%         of rows equal to the number of components. For the most general
%         case, this is equal to z.
%   y   - Initial guess for vapor mole fractions as a matrix with number
%         of rows equal to the number of components. For the most general
%         case, this is equal to z.
%
% OPTIONAL PARAMETERS:
%
%  'tol_equil'   - Tolerance for equilibrium constants in stability test.
%                  Uses the natural logarithm of K_i squared as a norm.
%
%  'tol_trivial' - Tolerance for the trivial test
%
% RETURNS:
%   stable - Indicator. If the entry for cell i is true, then the mixture
%            is stable (single-phase) for the p, T regime in that cell.
%   x      - Estimates of liquid mole fractions. Equal to z for stable
%            cells.
%   y      - Estimates of vapor mole fractions. Equal to z for stable
%            cells.
%
% SEE ALSO:
%   `EquationOfStateModel`

%{
Copyright 2009-2017 SINTEF Digital, Mathematics & Cybernetics.

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
    opt = struct('tol_equil', 1e-10, 'tol_trivial', 1e-5);
    assert(all(p > 0 & isfinite(p)), 'Positive, finite pressures required for phase stability test.');
    if numel(varargin)
        opt = merge_options(opt, varargin{:});
    end
    if nargin < 5
        x = z;
        y = z;
    end
    S_tol = opt.tol_trivial;
    nc = numel(p);
    active = true(nc, 1);
    
    z = ensureMinimumFraction(z);
    [y, S_V, isTrivialV] = checkStability(eos, z, y, p, T, true, active, opt);
    V_stable = S_V <= 1 + S_tol | isTrivialV;
    
    [x, S_L, isTrivialL] = checkStability(eos, z, x, p, T, false, V_stable, opt);
    L_stable = S_L <= 1 + S_tol | isTrivialL;

    stable = (L_stable & V_stable);
    x(stable, :) = z(stable, :);
    y(stable, :) = z(stable, :);
    
    bad = ~isfinite(S_L);
    if any(bad)
        for i = 1:numel(x)
            x(bad, :) = z(bad, :);
        end
    end
    bad = ~isfinite(S_V);
    if any(bad)
        for i = 1:numel(y)
            y(bad, :) = z(bad, :);
        end
    end
end

function [xy, S, trivialSolution, K] = checkStability(eos, z, xy, p, T, insidePhaseIsVapor, active, opt)
    tol_trivial = opt.tol_trivial;
    tol_equil = opt.tol_equil;
    nc = numel(p);

    [A_ij, Bi] = eos.getMixingParameters(p, T, eos.fluid.acentricFactors, false);
    
    f_z = getFugacity(eos, A_ij, Bi, z, p, T, insidePhaseIsVapor);
    K = estimateEquilibriumWilson(eos, p, T);
    % xy is either x or y, depending on context for phase test
    S = zeros(nc, 1);
    trivialSolution = false(nc, 1);
    if ~any(active)
        return
    end
    for stepNo = 1:20000
        zi = z(active, :);
        ki = K(active, :);
        A_ij_loc = cellfun(@(x) x(active, :), A_ij, 'UniformOutput', false);
        Bi_loc = Bi(active, :);
        
        if insidePhaseIsVapor
            xy_loc = ki.*zi;
        else
            xy_loc = zi./ki;
        end
        S_loc = sum(xy_loc, 2);
        xy_loc = bsxfun(@rdivide, xy_loc, S_loc);
        
        f_xy = getFugacity(eos, A_ij_loc, Bi_loc, xy_loc, p(active), T(active), ~insidePhaseIsVapor);

        f_zi = f_z(active, :);
        if insidePhaseIsVapor
            R = bsxfun(@rdivide, f_zi./f_xy, S_loc);
        else
            R = bsxfun(@times, f_xy./f_zi, S_loc);
        end
        K(active, :) = K(active, :).*R;
        R_norm = sum((R-1).^2, 2);
        K_norm = sum(log(K(active, :)).^2, 2);

        trivial = K_norm < tol_trivial;
        converged = R_norm < tol_equil;
        done = trivial | converged;
        S(active) = S_loc;
        trivialSolution(active) = trivial;
        xy(active, :) = xy_loc;
        active(active) = ~done;
        if all(done)
            dispif(mrstVerbose() > 1, 'Stability done in %d iterations\n', stepNo);
            break
        end
    end
    if ~all(done)
        warning('Stability test did not converge');
    end
end

function f = getFugacity(model, A_ij, Bi, xy, p, T, isLiquid)
    [Si, A, B] = model.getPhaseMixCoefficients(xy, A_ij, Bi);
    if isLiquid
        Z = model.computeLiquidZ(A, B);
    else
        Z = model.computeVaporZ(A, B);
    end
    f = model.computeFugacity(p, xy, Z, A, B, Si, Bi);
end