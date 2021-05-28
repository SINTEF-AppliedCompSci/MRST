function [stable, x, y] = phaseStabilityTest(eos, z, p, T, K, varargin)
% Perform a phase stability test for a mixture
%
% SYNOPSIS:
%   [stable, x, y] = phaseStabilityTest(eos, z, p, T, K)
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
%   K   - Initial guess for equilibrium constants. Will be estimated with
%         Wilson's correlation if K = [].
%
% OPTIONAL PARAMETERS:
%
%  'tol_equil'   - Tolerance for equilibrium constants in stability test.
%                  Uses the natural logarithm of K_i squared as a norm.
%
%  'tol_trivial' - Tolerance for the trivial test
%
%  'solve_both'  - Perform a stability test for a second liquid-like phase
%                  for cells already determined to be unstable by the
%                  second vapor-like phase. This has slightly higher
%                  computational cost since the stability must be checked
%                  twice, but the phase fraction estimates are more
%                  accurate. Enabled by default if x, y are requested as
%                  outputs.
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
    opt = struct('tol_equil',       1e-10, ...
                 'alpha',           [], ...
                 'tol_trivial',     1e-5, ...
                 'maxIterations',   20000, ...
                 'solve_both',      nargout > 1);
    assert(all(p > 0 & isfinite(p)), 'Positive, finite pressures required for phase stability test.');
    if numel(varargin)
        opt = merge_options(opt, varargin{:});
    end
    mc = eos.minimumComposition;
    if isempty(K)
        K = estimateEquilibriumWilson(eos, p, T);
    end

    S_tol = opt.tol_trivial;
    nc = numel(p);
    active = true(nc, 1);
    
    z = ensureMinimumFraction(z, mc);
    [y, S_V, isTrivialV] = checkStability(eos, z, K, p, T, true, active, opt);
    V_stable = S_V <= 1 + S_tol | isTrivialV;
    
    if ~opt.solve_both
        active = V_stable;
    end
    [x, S_L, isTrivialL] = checkStability(eos, z, K, p, T, false, active, opt);
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

function [xy, S, trivialSolution, K] = checkStability(eos, z, K, p, T, insidePhaseIsVapor, active, opt)
    tol_trivial = opt.tol_trivial;
    tol_equil = opt.tol_equil;
    nc = numel(p);
    acf = eos.CompositionalMixture.acentricFactors;
    [A_ij, Bi] = eos.getMixingParameters(p, T, acf, false);
    
    f_z = getFugacity(eos, A_ij, Bi, z, p, insidePhaseIsVapor);
    % xy is either x or y, depending on context for phase test
    S = zeros(nc, 1);
    trivialSolution = false(nc, 1);
    if ~any(active)
        return
    end
    % Init local parameters for eos
    A_ij_loc = cellfun(@(x) x(active, :), A_ij, 'UniformOutput', false);
    p_loc = p(active);
    z_loc = z(active, :);
    Ki = K(active, :);
    Bi_loc = Bi(active, :);
    f_zi = f_z(active, :);
    if insidePhaseIsVapor
        xy = K.*z;
    else
        xy = z./K;
    end
    for stepNo = 1:opt.maxIterations
        if insidePhaseIsVapor
            xy_loc = Ki.*z_loc;
        else
            xy_loc = z_loc./Ki;
        end
        S_loc = sum(xy_loc, 2);
        xy_loc = bsxfun(@rdivide, xy_loc, S_loc);
        
        f_xy = getFugacity(eos, A_ij_loc, Bi_loc, xy_loc, p_loc, ~insidePhaseIsVapor);
        if insidePhaseIsVapor
            % f_xy is vapor fugacity
            R = bsxfun(@rdivide, f_zi./f_xy, S_loc);
        else
            % f_xy is liquid fugacity
            R = bsxfun(@times, f_xy./f_zi, S_loc);
        end
        Ki = Ki.*R;
        R_norm = sum((R-1).^2, 2);
        K_norm = sum(log(Ki).^2, 2);
        
        trivial = K_norm < tol_trivial;
        converged = R_norm < tol_equil;
        done = trivial | converged;
        keep = ~done;
        % Store converged entries
        replace = active;
        replace(active) = done;
        K(replace, :) = Ki(done, :);
        S(replace) = S_loc(done, :);
        trivialSolution(replace) = trivial(done);
        xy(replace, :) = xy_loc(done, :);
        % Update active flag
        active(active) = ~done;
        
        if all(done)
            dispif(mrstVerbose() > 1, 'Stability done in %d iterations\n', stepNo);
            break
        end
        % Update local vectors for eos parameters
        A_ij_loc = cellfun(@(x) x(keep, :), A_ij_loc, 'UniformOutput', false);
        p_loc = p_loc(keep);
        z_loc = z_loc(keep, :);
        Ki = Ki(keep, :);
        Bi_loc = Bi_loc(keep, :);
        f_zi = f_zi(keep, :);
    end
    if ~all(done)
        warning('Stability test did not converge');
    end
end

function f = getFugacity(model, A_ij, Bi, xy, p, isLiquid)
    [Si, A, B] = model.getPhaseMixCoefficients(xy, A_ij, Bi);
    Z = model.computeCompressibilityZ(p, xy, A, B, Si, Bi, isLiquid);
    f = model.computeFugacity(p, xy, Z, A, B, Si, Bi);
end