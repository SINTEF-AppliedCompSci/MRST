function p_crit = computeCriticalPressure(relHumidity, ambientTemp, rho_w)
% Determine critical water pressure according to Simunek et al (2005)
% https://www.pc-progress.com/Downloads/Pgm_hydrus1D/HYDRUS1D-4.17.pdf 
%
% SYNOPSIS:
%   p_crit = criticalPressure(relHumidity, ambientTemp, rho_w)
%
% PARAMETERS:
%   relHumidity     - Scalar, relative humidity [%]
%   ambientTemp     - Scalar, ambient temperature [Celsius]
%   rho_w           - Scalar, water density [kg/m^3]
%
%  RETURNS:
%   p_cri           - Scalar, critical water pressure [Pa]
%

%{
Copyright 2018-2019, University of Bergen.

This file is part of the fv-unsat module.

fv-unsat is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

fv-unsat is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this file.  If not, see <http://www.gnu.org/licenses/>.
%} 

R = 8.314; % [J mol^{-1} K^{-1}] universal gas constant
temperature = ambientTemp + 273.15; % [K] absolute ambient temperature
M = 0.018015; % [kg mol^{-1}] molecular weight of water
p_crit = log(relHumidity/100) * R * temperature * rho_w / M;

end