function f = assignGRAVITY(f, api, reg)
%Define Fluid Densites from API and Specific Fluid Gravity
%
% SYNOPSIS:
%   fluid = assignGRAVITY(fluid, API, reg)
%
% PARAMETERS:
%   fluid - Partially formed ADI fluid object.
%
%   api   - Fluid gravity data comprising oil API degrees, specific gravity
%           of water relative to pure water, and specific gravity of gas
%           relative to air.
%
%   reg   - Region mapping structure.
%
% RETURNS:
%   fluid - Partitally formed ADI fluid object, now with added fluid
%           component mass density at surface conditions.
%
% NOTE:
%   The American Petroleum Institute defines the API gravity of oil as
%
%      API = (141.5 / SG) - 131.5
%
%   with 'SG' representing the specific gravity of oil relative to pure
%   water, measured at standard/surface conditions; 1 atmosphere and 60
%   degrees Fahrenheit.  API = 10 corresponds to an oil component whose
%   mass density is equal to that of pure water at surface conditions.
%
% SEE ALSO:
%   `initDeckADIFluid`.

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

   rho0 = reference_density();
   rho  = struct('O', oil_density(rho0.wat, api(:,1)), ...
                 'W', rho0.wat .* api(:,2),            ...
                 'G', rho0.air .* api(:,3));

   if numel(reg.PVTINX) == 1
      f.rhoOS = rho.O(1);
      f.rhoWS = rho.W(1);
      f.rhoGS = rho.G(1);
   else
      f.rhoOS = rho.O(reg.PVTNUM);
      f.rhoWS = rho.W(reg.PVTNUM);
      f.rhoGS = rho.G(reg.PVTNUM);
   end
end

%--------------------------------------------------------------------------

function rho = reference_density()
% At standard conditions--1 Atm & 60 Fahrenheit (15.56 Celsius); E100.
   rho = struct('wat', 1000*kilogram/meter^3, ...
                'air', 1.22*kilogram/meter^3);
end

%--------------------------------------------------------------------------

function rho = oil_density(rhoW0, api)
% Degrees API defined by American Petroleum Institute as
%
%    API = (141.5 / SG) - 131.5
%
% with SG = specific gravity of oil relative to pure water.

   rho = rhoW0 .* (141.5 ./ (api + 131.5));
end
