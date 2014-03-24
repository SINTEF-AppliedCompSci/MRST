function bo_fluid = fluid2BlackOilT(bo_fluid,fluids,p_ref,T_ref)
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

