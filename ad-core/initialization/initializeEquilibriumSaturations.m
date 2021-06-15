function s = initializeEquilibriumSaturations(model, region, pressures)
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

    nph = size(pressures, 2);
    refIx = region.reference_index;
    
    s = zeros(size(pressures));
        
    p_ref = pressures(:, refIx);
    for i = 1:nph
        if i == refIx
            continue
        end
        p = pressures(:, i);
        s(:, i) = solveSaturations(p, p_ref, ...
                         region.pc_functions{i}, region.pc_sign(i), ...
                         region.pc_scale(:, i), ...
                         region.s_min(:, i), region.s_max(:, i));
    end
    s(:, refIx) = 1 - sum(s, 2);
    
    bad = s(:, refIx) < 0;
    
    if any(bad)
        % Skip the middle phase and take the total capillary pressure to
        % find saturations
        s(bad, :) = 0;
        first = 1;
        last = nph;
        assert(refIx ~= first & refIx ~= last);
        % This may have the wrong sign.
        pc = @(S) (region.pc_sign(first)*region.pc_functions{first}(S) - ...
                   region.pc_sign(last)*region.pc_functions{last}(S));

        p_ref = pressures(bad, first);
        p = pressures(bad, last);
        
        s(bad, last) = solveSaturations(p, p_ref, ...
                         pc, 1, region.pc_scale(bad, 1),...
                         region.s_min(:, last), region.s_max(:, last));
        s(bad, first)= 1 - s(bad, last);
    end
end

function s = solveSaturations(p, p_ref, pc_fn, pc_sign, pc_scale, s_min, s_max)
    s = zeros(size(p));
    sat = (0:0.01:1)';
    pc = pc_sign*pc_fn(sat);

    dp =  p - p_ref;
    dp = dp./pc_scale;
    toMax = dp > max(pc);
    toMin = dp <= min(pc);
    if size(s_min, 1) == 1
        s(toMin) = s_min;
    else
        s(toMin) = s_min(toMin);
    end
    if size(s_max, 1) == 1
        s(toMax) = s_max;
    else
        s(toMax) = s_max(toMax);
    end
    middle = ~(toMin | toMax);
    if any(middle)
        s_inv = invertCapillary(dp(middle), pc_fn, pc_sign);
        fm = find(middle);
        bad = isnan(s_inv);
        if any(bad)
            if size(s_max, 1) == 1
                s_inv(bad) = s_max;
            else
                s_inv(bad) = s_max(fm(bad));
            end
        end
        if numel(s_min) == 1
            s_min = repmat(s_min, size(s_inv));
        else
            s_min = s_min(middle);
        end
        s(middle) = max(s_inv, s_min);
    end
end

function S = invertCapillary(dp, pc_fn, sgn)
    s = (0:0.0001:1)';
    pc = sgn*pc_fn(s);
    if pc(1) == pc(end)
        % There's no capillary pressure to speak of. Return NaN and fix
        % this on the outside.
        S = nan(size(dp));
        return;
    end
    
    [pc, ix] = unique(pc, 'last');
    s = s(ix);
    
    T = griddedInterpolant(pc, s, 'linear', 'nearest');
    S = T(dp);
end
