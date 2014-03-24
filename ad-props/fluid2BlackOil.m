function fluid = fluid2BlackOil(bo_fluid,fluids,p_ref,T_ref)
  %bo_fluid = fluid2BlackOil(bo_fluid,fluids,p_ref,T_ref)
  % relperm needs to be defined on bo_fluids
    names = {'W', 'O', 'G'};
   for i = 1:numel(names)
       n = names{i};       
       fluid.(['rho', n, 'S']) =  fluids{i}.density(p_ref,T_ref);
       fluid.(['b', n]) =@(p) fluids{i}.density(p,extend(T_ref,p))./fluid.(['rho', n, 'S']);
       fluid.(['B', n]) =@(p) 1./fluid.(['b', n])(p);
       try
            [pa,Ta] = initVariablesADI(p_ref,T_ref);
            fluids{i}.viscosity(pa,Ta);
            fluid.(['mu', n]) = @(p) fluids{i}.viscosity(p,extend(T_ref,p));
       catch
            warning('Derivatives of viscosity not present: use constant')
            fluid.(['mu', n]) = @(p) fluids{i}.viscosity(p_ref,T_ref);
       end
       fluid.(['h', n]) =@(p,T) fluids{i}.enthalpy(p,T);
       fluid.(['u', n]) =@(p,T) fluids{i}.enthalpy(p,T)-p./fluids{i}.density(p,T);
       fluid.(['kr',n]) = bo_fluid.(['kr',n]);
       fluid.(['krO', n]) = bo_fluid.(['kr',n]);
   end
   fluid.rsSat = @(varargin) varargin{1}*0;
   fluid.relPerm = bo_fluid.relPerm;
end
function T=extend(T,p)
    T=repmat(T,numel(double(p)),1);
end
