function K = estimateEquilibriumWilson(eos, p, T)
% Estimate equilibrium constant for a given pressure and temperature
    acf = eos.fluid.acentricFactors;
    % Estimate equilibrium constants using Wilson equation
    [Pr, Tr] = eos.getReducedPT(p, T, false);
    K = exp(5.37.*(bsxfun(@times, 1 + acf, 1 - 1./Tr)))./Pr;
    K(~isfinite(K)) = 1e3;
end