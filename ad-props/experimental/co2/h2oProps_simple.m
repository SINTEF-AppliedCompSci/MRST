function obj=h2oProps_simple()
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

warning('Incompressible water with pressure independent enthalpy')
obj.density =@(p,T) 1000+0*p;
obj.enthalpy =@(p,T)  4180*(T - 293.15)+0*p;
obj.viscosity = @(p,T)  1e-03+0*p;
obj.name='h2o simple';
end
