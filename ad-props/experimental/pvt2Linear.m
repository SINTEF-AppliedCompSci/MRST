function obj= pvt2Linear(fluid,p_ref,T_ref)
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

    [p,T]=initVariablesADI(p_ref,T_ref);
    
    rho=fluid.density(p,T);
    cp_rho=rho.jac{1}(1,1);ct_rho=rho.jac{2}(1,1);
    rho=double(rho);
    assert(cp_rho>=0);
    try
    mu=fluid.viscosity(p,T);
    cp_mu=mu.jac{1}(1,1);ct_mu=mu.jac{2}(1,1);
    mu=double(mu);
    catch
       warning('constant viscosity')
       mu=fluid.viscosity(double(p),double(T));
       cp_mu=0;ct_mu=0;
    end

    h=fluid.enthalpy(p,T);
    cp_h=full(h.jac{1}(1,1));ct_h=full(h.jac{2}(1,1));
    assert(ct_h>=0)
    h=double(h);
    if (abs((1 + ct_rho.*T_ref/rho)./rho - cp_h)>1e-3*abs(cp_h))
        warning('Fluid and enthalpy not consistent');
    end
    
   obj=linearFluid(cp_rho,ct_rho,ct_h,cp_mu,ct_mu,rho,h,mu,T_ref,p_ref);   

end
