function obj= pvt2Linear(fluid,p_ref,T_ref);
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