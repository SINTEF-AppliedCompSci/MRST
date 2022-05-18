function muB = viscosity_brine(p, t, c)
%Compute viscosity of a brine following Spivey et al., 2004
%
% SYNOPSIS:
%   mu_b = viscosity_brine(p,t,c)
% 
% PARAMETERS: 
%   p   - pressure. One value per cell. 
%   t   - temperature. One value per cell. 
%   c   - NaCl mass fraction. One value per cell.
% 
% RETURNS:
%   mu_b - brine viscosity. One value per cell 
% 
% SEE ALSO:
%   'density_pure_water', 'viscosity_pure_water', 'density_brine'

%{
Copyright 2009-2022 SINTEF ICT, Applied Mathematics.

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

    % Determine verbosity
    verbose = mrstVerbose();
    % Pressure
    pMinEOS = 0;       % Minimum pressure
    pMaxEOS = 200*1e6; % Maximum pressure
    % Trim pressure to EOS range
    p = trimToEOSRange(p, pMinEOS, pMaxEOS, 'pressure', verbose);
    p = p./1e6;  % conversion from Pa tp MPa 
    % Temperature
    t       = t-273.15; % Conversion from Kelvin to Celsius
    tMinEOS = 1e-3;     % Minimum temperature
    tMaxEOS = 275;      % Maximum temperature
    % Trim temperature to EOS range
    t = trimToEOSRange(t, tMinEOS, tMaxEOS, 'temperature', verbose);
    % NaCl molality
    mwNaCl        = 0.0584428;     % [kg/mol] Molar weigth of NaCl
    XNaCl         = c;             % Mass fraction of NaCl  
    XW            = 1 - XNaCl;     % Mass fraction of water    
    cNaCl         = XNaCl./XW;     % Concentration of NaCl per kg of water [kg/kg]   
    molNaCl       = cNaCl./mwNaCl; % Molality [mol/kg H2O]
    molNaClMinEOS = 0;             % Minimum molality
    molNaClMaxEOS = 5.7;           % Maximum molality
    % Trim NaCl molality to EOS range
    molNaCl = trimToEOSRange(molNaCl, molNaClMinEOS, molNaClMaxEOS, 'mol NaCl', verbose);
    % Coefficient for the scaling
    A1 = 0.0173; 
    A2 = 0.068;
    B1 = -1.0531; 
    B2 = 0.0273;
    % Function kestin_brine_viscosity takes MPa in input
    muKestin = kestin_brine_viscosity(p, 125,molNaCl);
    muK      = kestin_brine_viscosity(p,t,molNaCl);
    term1    = (A2.*(p./100)+A1).*(log(t./125)).^2;
    term2    = (B2.*(p./100)+B1).*(log(t./125));
    mu       = muKestin.*exp(term1 + term2);
    % Use Kestin values for t < 125
    useKestin = value(t) < 125;
    muB = mu.*(~useKestin) + muK.*(useKestin);
    muB = muB.*1e-6;

end

%--------------------------------------------------------------------------
%                             References
%--------------------------------------------------------------------------
% Spivey, J.P., Mccain, W.D., North, R. 2004. Estimating density, formation
% volume factor, compressibility, methane solubility, and viscosity for
% oilfield brines at temperatures from 0 to 275ºC, pressures to 200 Mpa,
% and salinities to 5.7 mol.kg-1, Journal of Canadian Petroleum Technology,
% 43 (7), 52-60.