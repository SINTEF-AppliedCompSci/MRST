function [poro, C]  = poroFromPerm(perm, avg_poro, tol)

% Compute porosity from permeability using the inverse Cozeny-Karman
% relation: K = C φ^3 / (1-φ)^2

% Compute C in order to target the correct average
avg_perm = mean(perm(:));
C = avg_perm * (1 - avg_poro)^2 / avg_poro^3;


% define function to find root of, as well as its derivative
fac = perm / C;
f = @(phi) phi.^3 - fac .* (1-phi).^2;
df = @(phi) 3 * phi.^2 + 2 * fac .* (1-phi);

% initial guess for poro
poro = avg_poro * ones(size(perm)); % perm should be a vector here

% iterate until poro satisfies relationship
[res, d_res] = deal(f(poro), df(poro));

while max(abs(res)) > tol 
    d_phi = res ./ d_res;
    poro = poro - d_phi;
    
    [res, d_res] = deal(f(poro), df(poro));
end
