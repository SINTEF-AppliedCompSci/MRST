function obj= pvt2Linear(fluid,p_ref,T_ref);
    [p,T]=initVariablesADI(p_ref,T_ref);
    rho=fluid.density(p,t);
    cp_rho=rho.jac{1}(1,1);ct_rho=rho.jac{2}(1,1);
    rho=double(rho)
    rho=fluid.density(p,t);
    cp_rho=rho.jac{1}(1,1);ct_rho=rho.jac{2}(1,1);
    rho=double(rho)
    
   obj=linearFluid(cp_rho,ct_rho,ct_h,cp_mu,ct_mu,rho,mu,T_ref,p_ref)   

end