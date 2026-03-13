function tab_sol = HenrySetschenowH2BrineEos(T, mNaCl, p)
% HENRYSETSCHENOW: Calculates the mole fraction of dissolved H2 (xH2)
% in solution using the Henry–Setschenow equation with temperature and
% salinity dependence.
%
% INPUTS:
% T      - Temperature in Kelvin (K)
% mNaCl  - Salinity in mol/kg (molal concentration of NaCl)
% p      - Partial pressure of H2 in Pascal
%
% OUTPUT:
% xH2    - Mole fraction of H2 in solution
% REFRENCE
%
%   Ahmed, E., et al. (2024). Phase behavior and black-oil simulations of 
%   hydrogen storage in saline aquifers. Advances in Water Resources, 191, 104772.
%
%{
Copyright 2009-2025 SINTEF Digital, Mathematics & Cybernetics.

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
% Constants for the temperature-dependent Henry’s law calculation
kH0 = 1282.05.*atm*litre*mol^-1;  % [atm L mol^-1]
kH1 = 1.8018.*mol.*atm^-1;   % [mol atm^-1]
a = 6.9376;
b = -7.1384;
c = -12.6;
d = -2.7320;

% Constants for the salinity dependence
e = 0.0185;
f = -0.3019;

% Calculate reduced temperature
Tr = T / 647.096;  % Reference temperature (298.15 K)

% Calculate log(kH) based on temperature dependence
log_kH = log(kH0 ./ kH1) + a + b ./ Tr + c * log(Tr) + d .* Tr;
kH_T = exp(log_kH);  % Henry's constant kH(T)

% Calculate log(qH) based on salinity dependence
log_qH = e * mNaCl.^2 + f .* mNaCl;
qH = exp(log_qH);  % Salinity correction factor qH(mNaCl)

% Calculate the mole fraction of H2 in solution using Henry–Setschenow equation
xH2 = (qH .* p) ./ kH_T;
T = T(:);       % ensure column vector
p = p(:);       % ensure column vector
xH2 = xH2(:);   % ensure column vector

tab_sol = table();
tab_sol.("# temperature [°C]") = T - 273.15;
tab_sol.("pressure [Pa]") = p;
tab_sol.("x_H2") = xH2;


end
