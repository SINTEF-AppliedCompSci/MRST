function s = initializeEquilibriumSaturations(model, region, pressures)
    nph = size(pressures, 2);
    refIx = region.reference_index;
    
    s = zeros(size(pressures));
        
    p_ref = pressures(:, refIx);
    for i = 1:nph
        if i == refIx
            continue
        end
        
%         pc = region.pc_sign(i)*region.pc_functions{i}(sat);
        p = pressures(:, i);
        s(:, i) = solveSaturations(p, p_ref, ...
                         region.pc_functions{i}, region.pc_sign(i),...
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
                         pc, 1,...
                         region.s_min(:, last), region.s_max(:, last));
        s(bad, first)= 1 - s(bad, last);
    end
end

function s = solveSaturations(p, p_ref, pc_fn, pc_sign, s_min, s_max)
    s = zeros(size(p));
    
    sat = 0:0.01:1;
    pc = pc_sign*pc_fn(sat);
    dp =  p - p_ref;

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
        s(middle) = s_inv;
    end
end

function S = invertCapillary(dp, pc_fn, sgn)
    s = (0:0.0001:1)';
    pc = sgn*pc_fn(s);
    
    [pc, ix] = unique(pc, 'last');
    s = s(ix);
    
    T = griddedInterpolant(pc, s, 'linear', 'nearest');
    S = T(dp);
end
