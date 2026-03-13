function Z = calculateBrillBreggsZfactorHydrogen(T, P)
% calculateBrillBeggsZfactorHydrogen  Computes the hydrogen density using the Brill and Beggs correlation.
%
% SYNOPSIS:
%   rho = calculateBrillBeggsZfactorHydrogen(T, P)
%
% DESCRIPTION:
%   This function calculates the hydrogen density based on the Brill and Beggs
%   correlation, incorporating temperature (T) and pressure (P). The correlation
%   is derived from empirical coefficients specific to hydrogen gas properties.
%
% INPUTS:
%   T - Temperature in Kelvin (K)
%   P - Pressure in megapascals (MPa)
%
% OUTPUTS:
%   rho - Hydrogen density in kilograms per cubic meter (kg/m^3)
%
% REFERENCE:
%   Jafari Raad, Seyed Mostafa, Leonenko, Yuri, Hassanzadeh, Hassan, 2023.
%   Correlations for prediction of hydrogen gas viscosity and density for
%   production, transportation, storage, and utilization applications.
%   Int. J. Hydrog. Energy (ISSN: 0360-3199) 48 (89), 34930–34944.

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


Tc = 33.19; % Critical temperature in Kelvin
pc = 1.315e6; % Critical pressure in Pascal

% Calculate pseudo-reduced temperature and pressure
Tpr = T ./ Tc; % Pseudo-reduced temperature
ppr = P ./ pc; % Pseudo-reduced pressure

% Coefficients
c1 = 0.4492195; c2 = -23.84114; c3 = 0.0562519;
c4 = 1.01124; c5 = -0.0289739; c6 = -0.0011064;
c7 = -0.0385377; c8 = -42.98531; c9 = 0.0010446;
c10 = 0.004043; c11 = 0.0043762; c12 = -70.18049;
c13 = 0.016152; c14 = 0.0092372; c15 = 0.1085783;
c16 = 0.0076386; c17 = -9.09E-05;

% Calculate parameters A, B, C, and D
A = c1 .* (Tpr - c2).^(0.5) - c3 .* Tpr - c4;
B = (c5 - c6 .* Tpr).* ppr + (c7./(Tpr - c8) - c9 + c10 ./ (10.^(c11 .* (Tpr - c12)))).* ppr.^2;
C = c13 - c14 * log(Tpr);
D = 10.^(c15) * c16 .* Tpr + c17 .* Tpr.^2;

% Calculate Z using the provided formula
Z = A + (1 - A)./exp(B) + C.* ppr.^D;
end
