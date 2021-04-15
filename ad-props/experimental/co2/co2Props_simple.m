function obj = co2Props_simple()
% from dumux

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

warning('CO2 as ideal gas')
obj.density =@(p,T) g_density(T,p);
obj.enthalpy =@(p,T) g_enthalpy(T,p);
obj.viscosity = @(p,T) g_viscosity(T,p);
obj.name='co2 simple';
end

function v= l_enthalpy(T,p)
v= (T - 298.15)*5e3;
end

function v = g_enthalpy(T,p)
   v = 571.3e3 + (T - 298.15)*0.85e3;
end

function v=g_density(T,p)
mM=44e-3;%molarMass()
R = 8.314472;
%idal gass
 v = p.*mM./(R*T);
end
function v= g_viscosity(T,p)

Tc = 273.15 + 30.95;%criticalTemperature();
Vc = 93.9; %// critical specific volume [cm^3/mol]
omega = 0.239; %// accentric factor
mM=44e-3;%molarMass()
M = mM * 1e3; %// molar mas [g/mol]
dipole = 0.0; %// dipole moment [debye]

mu_r4 = 131.3 * dipole / sqrt(Vc * Tc);
mu_r4= power(mu_r4,3);
%mu_r4 *= mu_r4;
%mu_r4 *= mu_r4;

Fc = 1 - 0.2756*omega + 0.059035*mu_r4;
Tstar = 1.2593 * T/Tc;
Omega_v =...
    1.16145*power(Tstar, -0.14874) +...
    0.52487*exp(- 0.77320*Tstar) +...
    2.16178*exp(- 2.43787*Tstar);
mu = 40.785*Fc.*power(M*T,0.5)./(power(Vc, 2./3)*Omega_v);

%// convertion from micro poise to Pa s
v=mu/1e6 / 10;
end
