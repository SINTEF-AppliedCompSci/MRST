function muW = viscosity_pure_water(p, t)
%Compute viscosity of pure water following Spivey et al., 2004
%
% SYNOPSIS:
%   mu_w = viscosity_pure_water(p,t)
% 
% PARAMETERS: 
%   p   - pressure. One value per cell. 
%   t   - temperature. One value per cell. 
% 
% RETURNS:
%   mu_w - pure water viscosity. One value per cell 
% 
% SEE ALSO:
%   'density_brine', 'density_pure_water', 'viscosity_brine'

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
    
    % Call brine viscosity function with zero NaCl mass fraction
    muW = viscosity_brine(p, t, 0);
    
end

%--------------------------------------------------------------------------
%                             References
%--------------------------------------------------------------------------
% Spivey, J.P., Mccain, W.D., North, R. 2004. Estimating density, formation
% volume factor, compressibility, methane solubility, and viscosity for
% oilfield brines at temperatures from 0 to 275ºC, pressures to 200 Mpa,
% and salinities to 5.7 mol.kg-1, Journal of Canadian Petroleum Technology,
% 43 (7), 52-60.

