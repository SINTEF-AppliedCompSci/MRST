function bo_fluid = fluid2BlackOilT(bo_fluid,fluids,p_ref,T_ref)
%Undocumented Utility Function

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

  %bo_fluid = fluid2BlackOil(bo_fluid,fluids,p_ref,T_ref)
  % relperm needs to be defined on bo_fluids
    names = {'W', 'O', 'G'};
   for i = 1:numel(names)
       n = names{i};       
       fluid.(['rho', n, 'S']) =@(p,T) fluids{i}.density(p_ref,T_ref);
       fluid.(['b', n]) =@(p,T) fluids{i}.density(p,T)./fluid.(['rho', n, 'S']);
       fluid.(['B', n]) = 1./fluid.(['b', n]);
       fluid.(['mu', n]) = @(p,T) fluids{i}.viscosity(p,T);
       fluid.(['h', n]) =@(p,T) fluids{i}.enthalpy(p,T);
   end
   fluid.rsSat = @(varargin) varargin{1}*0;   
end

