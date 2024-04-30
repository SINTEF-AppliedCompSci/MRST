function [poro, C]  = poroFromPerm(perm, avg_poro, tol)

% Compute porosity from permeability using the inverse Cozeny-Karman
% relation: K = C φ^3 / (1-φ)^2

%{
Copyright 2009-2024 SINTEF Digital, Mathematics & Cybernetics.

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
