function [stable, x, y, L_stable, V_stable] = PhaseStabilityTest(eos, z, p, T)
    z = ensureMinimumFraction(z);

    [y, S_V, isTrivialV] = checkStability(eos, z, p, T, true);
    [x, S_L, isTrivialL] = checkStability(eos, z, p, T, false);
    
    L_stable = S_L <= 1 | isTrivialL;
    V_stable = S_V <= 1 | isTrivialV;
    stable = (L_stable & V_stable);
    
    bad = ~isfinite(S_L);
    if any(bad)
        for i = 1:numel(x);
            x{i}(bad) = z{i}(bad);
        end
    end
    bad = ~isfinite(S_V);
    if any(bad)
        for i = 1:numel(y);
            y{i}(bad) = z{i}(bad);
        end
    end
end

function [xy, S, trivialSolution] = checkStability(eos, z, p, T, insidePhaseIsVapor)
    ncomp = numel(z);
    f_z = getFugacity(eos, z, p, T, insidePhaseIsVapor);
    K = EstimateEquilibriumWilson(eos, p, T);
    % xy is either x or y, depending on context for phase test
    xy_loc = cell(1, ncomp);
    xy = cell(1, ncomp);
    nc = numel(p);
    [xy{:}] = deal(zeros(nc, 1));
    active = true(nc, 1);
    S = zeros(nc, 1);
    trivialSolution = false(nc, 1);
    for stepNo = 1:100
        S_loc = 0;
        for i = 1:ncomp
            zi = z{i}(active);
            ki = K{i}(active);
            if insidePhaseIsVapor
                xy_loc{i} = ki.*zi;
            else
                xy_loc{i} = zi./ki;
            end
            S_loc = S_loc + xy_loc{i};
        end
        for i = 1:ncomp
            xy_loc{i} = xy_loc{i}./S_loc;
        end
        f_xy = getFugacity(eos, xy_loc, p(active), T(active), ~insidePhaseIsVapor);
        
        R = cell(ncomp, 1);
        R_norm = 0;
        K_norm = 0;
        for i = 1:ncomp
            f_zi = f_z{i}(active);
            if insidePhaseIsVapor
                R{i} = (f_zi./f_xy{i}).*(1./S_loc);
            else
                R{i} = (f_xy{i}./f_zi).*S_loc;
            end
            K{i}(active) = K{i}(active).*R{i};
            R_norm = R_norm + (R{i} - 1).^2;
            K_norm = K_norm + (log(K{i}(active))).^2;
        end
        
        trivial = K_norm < 1e-5;
        converged = R_norm < 1e-10;
        done = trivial | converged;
        S(active) = S_loc;
        trivialSolution(active) = trivial;
        for i = 1:ncomp
            xy{i}(active) = xy_loc{i};
        end
        active(active) = ~done;
        if all(done)
            break
        end
        
    end
end

function K = EstimateEquilibriumWilson(eos, p, T)
    ncomp = eos.fluid.getNumberOfComponents();
    acf = eos.fluid.acentricFactors;
    % Estimate equilibrium constants using Wilson equation
    K = cell(1, ncomp);
    for i = 1:ncomp
        Tr = T./eos.fluid.Tcrit(i);
        Pr = p./eos.fluid.Pcrit(i);
        k = exp(5.37.*(1 + acf(i)).*(1 - 1./Tr))./Pr;
%         k(~isfinite(k)) = 1000;
        K{i} = k;
    end
    
end

function f = getFugacity(model, xy, p, T, isLiquid)
    [A_ij, Bi] = model.getMixingParameters(p, T, model.fluid.acentricFactors);
    [Si, A, B] = model.getPhaseMixCoefficients(xy, A_ij, Bi);
    if isLiquid
        Z = model.computeLiquidZ(A, B);
    else
        Z = model.computeVaporZ(A, B);
    end
    f = model.computeFugacity(p, xy, Z, A, B, Si, Bi);
end