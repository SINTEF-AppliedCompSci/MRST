function [foptval, uopt, history, uu_opt, extra] = ...
       optimize_mech(u, G, bcfun, cfun, loadfun, obj_fun)
%
%
% SYNOPSIS:
%   function [uopt, uu_opt, extra] = optimize_mech(u, G, bcfun, cfun, loadfun, obj_fun)
%
% DESCRIPTION:
%
% PARAMETERS:
%   u       - 
%   G       - 
%   bcfun   - 
%   cfun    - 
%   loadfun - 
%   obj_fun - 
%
% RETURNS:
%   foptval -
%   uopt    - 
%   history - 
%   uu_opt  - 
%   extra   - 
%
% EXAMPLE:
%
% SEE ALSO:
%   VEM_linelast_AD

   mrstModule add optimization;
   
   G = createAugmentedGrid(computeGeometry(G));
   
   funwrap = @(u) fun_wrapper(u, G, bcfun, cfun, loadfun, obj_fun);
   
   [foptval, uopt, history] = unitBoxBFGS(u, funwrap);

   %% compute additional information if requested
   if nargout > 3
      C = cfun(uopt);
      bc = bcfun(uopt);
      load = loadfun(uopt);

      if nargout == 4
         uu_opt = VEM_linElast_AD(G, C, bc, load);
      else
         [uu_opt, extra] = VEM_linElast_AD(G, C, bc, load);
      end
   end
end

function [val, grad] = fun_wrapper(u, G, bcfun, cfun, loadfun, obj_fun)

   u = initVariablesADI(u);
   
   bc = bcfun(u);
   C = cfun(u);
   load = loadfun(u);
   [dd, extra] = VEM_linElast_AD(G, C, bc, load);

   dofs = ~extra.disc.isdirdofs; %% exclude dirichlet nodes

   dd = dd';
   dd = dd(dofs);
   
   [val, oval_du, oval_dd] = obj_fun(u, dd(:));
   
   %% use adjoint to compute gradient
   lambda = -extra.A \ oval_dd; % A symmetric, so no transpose necessary
   
   dAdu_dd = 0; % @@ will change when including stiffness params. dependence on u
   dbdu = extra.rhs.jac{1};
   dsys_du = dAdu_dd - dbdu;
   
   grad = oval_du + lambda' * dsys_du;
   %grad = oval_du - lambda' * dsys_du;
      
end


   
   %function [v, u, history] = unitBoxBFGS(u0, f, varargin)   
%function [uu, extra] = VEM_linElast_AD(G, C, el_bc, load, varargin)
   
