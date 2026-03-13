function [stable, x, y, K] = eCPAphaseStabilityTest(eos, z, p, T, K, varargin)
% Perform a phase stability test for a mixture

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
Copyright 2009-2022 SINTEF Digital, Mathematics & Cybernetics.

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
        K = eCPAestimateEquilibriumWilson(eos, p, T);
    end

    S_tol = opt.tol_trivial;
    nc = numel(p);
    active = true(nc, 1);
    
    z = ensureMinimumFraction(z, mc);
    [y, S_V, isTrivialV, K_V] = checkStability(eos, z, K, p, T, true, active, opt);
    V_stable = S_V <= 1 + S_tol | isTrivialV;
    
    mu0 = eos.ECPACompositionalMixture.mu0;
    mu0 = mu0(~isnan(mu0));
    if ~opt.solve_both || ~isempty(mu0)
        active = V_stable;
    end
    [x, S_L, isTrivialL, K_L] = checkStability(eos, z, K, p, T, false, active, opt);
    L_stable = S_L <= 1 + S_tol | isTrivialL;
    
    K(~V_stable, :) = K_V(~V_stable, :);
    K(~L_stable, :) = K_L(~L_stable, :);
    
   
    stable = (L_stable & V_stable);
    x(stable, :) = z(stable, :);
    y(stable, :) = z(stable, :);
    x = x ./sum(x, 2);
    y = y ./sum(y, 2);

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
    nmole = eos.getNumberOfMolecules();
    
    % xy is either x or y, depending on context for phase test
    S = zeros(nc, 1);
    trivialSolution = false(nc, 1);
    if ~any(active)
        xy = z;
        return
    end
    
    K(:,nmole+1:end) = min(K(:,nmole+1:end), 1e-5);
    if insidePhaseIsVapor
        xy = K.*z;
    else
        z(:,nmole+1:end) = 1e-10;
        xy = z./K;
    end
    
    [A_ij, Ai, Tr] = eos.getMixingParameters(T, false);
    f_z = getFugacity(eos, T, z, p, A_ij, Ai, Tr, insidePhaseIsVapor);
    % Init local parameters for eos
    A_ij_loc = cellfun(@(x) x(active, :), A_ij, 'UniformOutput', false);
    p_loc = p(active);
    T_loc = T(active);
    z_loc = z(active, :);
    Ai_loc = Ai(active, :);
    Tr_loc = Tr(active, :);
    Ki = K(active, :);
    f_zi = f_z(active, :);
    for stepNo = 1:opt.maxIterations
        if insidePhaseIsVapor
            xy_loc = Ki.*z_loc;
        else
            xy_loc = z_loc./Ki;
        end
        xy_loc(:,nmole+1:end) = min(max(xy_loc(:,nmole+1:end), 1e-8), 0.01*xy_loc(:,1));
        S_loc = sum(xy_loc, 2);
        xy_loc = bsxfun(@rdivide, xy_loc, S_loc);

        f_xy = getFugacity(eos, T_loc, xy_loc, p_loc, A_ij_loc, Ai_loc, Tr_loc, ~insidePhaseIsVapor);
        if insidePhaseIsVapor
            % f_xy is vapor fugacity
            R = bsxfun(@rdivide, f_zi./f_xy, S_loc);
        else
            % f_xy is liquid fugacity
            R = bsxfun(@times, f_xy./f_zi, S_loc);
        end
        Ki = Ki .* R;
        R_norm = sum((R(:, 1:nmole) - 1).^2, 2);
        K_norm = sum(log(Ki(:, 1:nmole)).^2, 2);
        
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
        f_zi = f_zi(keep, :);
        T_loc = T_loc(keep);
        Ai_loc = Ai_loc(keep, :);
    end
    K(:,nmole+1:end) = min(max(K(:,nmole+1:end), 1e-5), 1);
    if ~all(done)
        warning('Stability test did not converge');
    end
end

function [f, v] = getFugacity(model, T, xy, p, A_ij, Ai, Tr, isLiquid)
    [Si, A, B] = model.getPhaseMixCoefficients(xy, A_ij);
    [~, v, XA, XC] = model.computeCompressibilityZ(p, T, xy, A, B, Ai, Si, Tr, isLiquid);
    f = model.computeFugacity(T, p, xy, A, B, Ai, v, Si, XA, XC, Tr);
end