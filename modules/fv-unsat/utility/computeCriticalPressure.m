function p_crit = computeCriticalPressure(phys)
% Determine critical water pressure according to Simunek et. al. (2005)
% https://www.pc-progress.com/Downloads/Pgm_hydrus1D/HYDRUS1D-4.17.pdf 
%
% SYNOPSIS:
%   p_crit = criticalPressure(phys)
%
% PARAMETERS:
%   phys            - Structure, containing the physical parameters
%
%  RETURNS:
%   p_crit          - Scalar, critical water pressure in Pascal.
%

%{
Copyright 2018-2020, University of Bergen.

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
M = 0.018015; % [kg mol^{-1}] molecular weight of water
p_crit = log(phys.flow.relativeHumidity) * R * phys.flow.temperature ...
    * phys.flow.rho / M;

end