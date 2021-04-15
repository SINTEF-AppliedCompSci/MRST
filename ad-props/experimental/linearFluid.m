function obj = linearFluid(cp_rho,ct_rho,ct_h,cp_mu,ct_mu,rho,h, mu,T_ref,p_ref)
%linearFluid(cp_rho,ct_rho,ct_h,cp_mu,ct_mu,rho,h, mu,T_ref,p_ref)

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

    %{
    obj.density =@(p,T) rho.*exp(cp_rho.*(p-p_ref)-ct_rho.*(T-T_ref));
    obj.viscosity =@(p,T) mu.*exp(-cp_mu.*(p-p_ref)+ct_mu.*(T-T_ref));
    %}
    obj.density =@(p,T) rho+cp_rho.*(p-p_ref)-ct_rho.*(T-T_ref);
    obj.viscosity =@(p,T) mu-cp_mu.*(p-p_ref)+ct_mu.*(T-T_ref);
    obj.enthalpy =@(p,T) enthalpy(p,T,p_ref,T_ref,ct_h,ct_rho,rho)+h;
    obj.name ='Linear fluid';
    %obj.energy  =@(p,T)  energy(p,T);                       
end
%{
function v=energy(p,T,obj,c_p,cp_rho,ct_rho)
    rho = obj.density(p,T);
    v= (c_p - ct_rho.*p./rho).*(T-T_ref)...
       + (cp_rho*p - ct_rho*T)./rho*(p-p_ref)
end
%}
function v=enthalpy(p,T,p_ref,T_ref,ct_h,ct_rho,rho)
    v= ct_h*(T-T_ref)...
       + ((1 + ct_rho.*T_ref/rho)./rho).*(p-p_ref);
end
