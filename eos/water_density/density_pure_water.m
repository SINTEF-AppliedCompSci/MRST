function rhoW = density_pure_water(p, t)
%Compute density of pure water following Spivey et al., 2004
%
% SYNOPSIS:
%   rho_b = density_pure_water(p,t)
% 
% PARAMETERS: 
%   p   - pressure. One value per cell. 
%   t   - temperature. One value per cell. 
% 
% RETURNS:
%   rho_w - pure water density. One value per cell 
% 
% SEE ALSO:
%   'density_brine', 'viscosity_pure_water', 'viscosity_brine'

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
    t  = t-273.15;   % Conversion from Kelvin to deg Celsius
    tMinEOS = 1e-3;  % Minimum temperature
    tMaxEOS = 275;   % Maximum temperaure
    % Trim temperature to EOS range
    t = trimToEOSRange(t, tMinEOS, tMaxEOS, 'temperature', verbose);
    % Get coefficients for pure water
    [rhoW0, Ew, Fw] = coefficients_pure_water( t );
    % Compute isothermal compressibility
    Iw  = isothermal_compressibility(p,p0,Ew,Fw);
    Iw0 = isothermal_compressibility(p0,p0,Ew,Fw);
    % Compute density of pure water
    rhoW = rhoW0 .* exp(Iw - Iw0); % eq. (8) Spivey et al., 2004
    rhoW = rhoW.*1000;
    
end

%--------------------------------------------------------------------------
%                             Reference
%--------------------------------------------------------------------------
% Spivey, J.P., Mccain, W.D., North, R. 2004. Estimating density, formation
% volume factor, compressibility, methane solubility, and viscosity for
% oilfield brines at temperatures from 0 to 275ºC, pressures to 200 Mpa,
% and salinities to 5.7 mol.kg-1, Journal of Canadian Petroleum Technology,
% 43 (7), 52-60.