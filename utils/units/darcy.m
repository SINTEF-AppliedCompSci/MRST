function d = darcy()
%Compute numerical value, in units of m^2, of the Darcy constant.
%
% SYNOPSIS:
%   d = darcy()
%
% PARAMETERS:
%   None.
%
% RETURNS:
%   d - Numerical value, in units of m^2, of the Darcy constant.
%
% NOTE:
%   A porous medium with a permeability of 1 darcy permits a flow (flux) of
%   1 cm³/s of a fluid with viscosity 1 cP (1 mPa·s) under a pressure
%   gradient of 1 atm/cm acting across an area of 1 cm².
%
% SEE ALSO:
%   `gravity`.

%{
Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.

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


   mu  = centi()*poise(); % Fluid viscosity, 1cP = 1mPa·s
   cm  = centi()*meter(); % [m]
   sec = second();        % [s]

   p_grad = atm() / cm;   % Pressure gradient        [Pa/m]
   area   = cm ^ 2;       % Active area              [m^2]
   flow   = cm ^ 3 / sec; % Flow rate                [m^3/s]
   vel    = flow / area;  % Fluid velocity           [m/s]

   d = vel * mu / p_grad; % == 1e-7 [m^2] / 101325
                          % == 9.869232667160130e-13 [m^2]
end
