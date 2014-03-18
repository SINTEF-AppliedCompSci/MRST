function obj = linearFluid(cp_rho,ct_rho,ct_h,cp_mu,ct_mu,rho,mu,T_ref,p_ref)   
    obj.density =@(p,T) rho.*exp(cp_rho.*(p-p_ref)-ct_rho.*(T-T_ref));
    obj.viscosity =@(p,T) mu.*exp(-cp_mu.*(p-p_ref)+ct_mu.*(T-T_ref));
    obj.enthalpy =@(p,T) enthalpy(p,T,p,T,obj,ct_h,ct_rho);
    obj.energy  =@(p,T)  energy(p,T);                       
end
function v=energy(p,T,obj,c_p,cp_rho,ct_rho)
    rho = obj.density(p,T);
    v= (c_p - ct_rho.*p./rho).*(T-T_ref)...
       + (cp_rho*p - ct_rho*T)./rho*(p-p_ref);
end
function v=enthalpy(p,T,obj,ct_h,ct_rho)
    rho = obj.density(p,T);
    v= ct_h*(T-T_ref)...
       + ((1 - ct_rho.*T)./rho).*(p-p_ref);
end