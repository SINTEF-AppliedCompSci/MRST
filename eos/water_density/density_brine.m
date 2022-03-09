function rho_b = density_brine(p, t, c)
%Compute density of a brine following Spivey et al., 2004
%
% SYNOPSIS:
%   rho_b = density_brine(p,t,c)
% 
% PARAMETERS: 
%   p   - pressure. One value per cell. 
%   t   - temperature. One value per cell. 
%   c   - NaCl mass fraction. One value per cell.
% 
% RETURNS:
%   rho_b - brine density. One value per cell 
% 
% SEE ALSO:
%   'density_pure_water', 'viscosity_pure_water', 'viscosity_brine'

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
    pMinEOS = 0;     % Maximum pressure
    pMaxEOS = 200e6; % Maximum pressure
    p0      = 70e6;  % [Pa] Reference pressure
    % Trim pressure to EOS range
    p = trimToEOSRange(p, pMinEOS, pMaxEOS, 'pressure', verbose);
    % Temperature
    t  = t-273.15; % Conversion from Kelvin to deg Celsius
    tMinEOS = 0;   % Minimum temperature
    tMaxEOS = 275; % Maximum temperaure
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

    % Get coefficients for pure water to compute density brine (rho_b0) and
    % coefficients for the brine (Eb, Fb)
    [rho_w0, Ew, Fw] = coefficients_pure_water( t );

    % Compute the density of the brine at the reference pressure
    rho_b0 = density_brine_refP(rho_w0,t,molNaCl);

    % Compute the coefficients Eb and Fb for the brine
    [Eb,Fb] = coefficients_brine(Ew,Fw,molNaCl,t);

    % Compute isothermal compressibility
     Ib  = isothermal_compressibility(p,p0,Eb,Fb);
     Ib0 = isothermal_compressibility(p0,p0,Eb,Fb);

    % Compute density
     rho_b = rho_b0 .* exp(Ib - Ib0); % Eq.(15) in Spivey et al., 2004
     rho_b = rho_b.*1000;

end

%--------------------------------------------------------------------------
%                             References
%--------------------------------------------------------------------------
% Spivey, J.P., Mccain, W.D., North, R. 2004. Estimating density, formation
% volume factor, compressibility, methane solubility, and viscosity for
% oilfield brines at temperatures from 0 to 275ºC, pressures to 200 Mpa,
% and salinities to 5.7 mol.kg-1, Journal of Canadian Petroleum Technology,
% 43 (7), 52-60.