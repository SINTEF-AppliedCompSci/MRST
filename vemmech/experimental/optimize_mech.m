function [foptval, uopt, history, uu_opt, extra] = ...
       optimize_mech(u, G, bcfun, cfun, loadfun, obj_fun, varargin)
%
%
% SYNOPSIS:
%   function [uopt, uu_opt, extra] = optimize_mech(u, G, bcfun, cfun, loadfun, obj_fun)
%
% DESCRIPTION:
%
% PARAMETERS:
%   u       - vector of control parameters (initial guess)
%   G       - simulation grid
%   bcfun   - function that generates an AD-struct with boundary conditions
%             from the control parameters
%   cfun    - function that generates an AD-tensor with material properties
%             from the control parameters (@@ AD currently unsupported)
%   loadfun - function that generates the cellwise load term (AD) from the
%             control parameters
%   obj_fun - objective function to minimize.  Should take the vector of
%             control variables as its first argument, the computed
%             displacements as its second argument, and return the following:
%             + val - the objective value
%             + du  - partial derivatives of the objective value wrt. control
%                     variables
%             + dd  - partial derivatives of the objectiv value wrt. 
%                     displacements
%
% RETURNS:
%   foptval - the optimal objective function value found
%   uopt    - the corresponding optimal control parameter values
%   history - a struct containing history about the optimization procedure,
%             including the intermediary choices of control variables and
%             objective function values during the search.
%   uu_opt  - displacements evaluated with the optimal control parameter values
%   extra   - the extra information returned by VEM_linElast_AD, including
%             system matrix, right-hand-side, discrete operators, etc.
%
% SEE ALSO:
%   VEM_linelast_AD

   mrstModule add optimization;
   
   opt.gradTol = 1e-3;
   opt.objChangeTol = 1e-5;
   opt.cyclical = []; % indices of cyclical control variables
   [opt, ~] = merge_options(opt, varargin{:});

   
   G = createAugmentedGrid(computeGeometry(G));
   
   funwrap = @(u) fun_wrapper(u, G, bcfun, cfun, loadfun, obj_fun);
                              
   
   [foptval, uopt, history] = unitBoxBFGS(u, funwrap, ...
                                          'gradTol', opt.gradTol, ...
                                          'objChangeTol', opt.objChangeTol);

   % handle cyclical variables that have hit against the imposed box boundary
   small = 1e-2;
   if ~isempty(opt.cyclical)
      
      ixs_1 = uopt(opt.cyclical) == 1;
      ixs_0 = uopt(opt.cyclical) == 0;
      
      while any([ixs_1(:) ;ixs_0(:)]) && ...
             (history.pg(end) > opt.gradTol || isnan(history.pg(end)))
         uopt
         % one or more cyclical variables 'stuck' against the imposed boundary.  Move them
         % to the other side and try again.
         fprintf('cyclical variable(s) stuck.  Trying again.\n');
         u = uopt;
         u(opt.cyclical(ixs_1)) = 0+small;
         u(opt.cyclical(ixs_0)) = 1-small;
         [foptval, uopt, history] = unitBoxBFGS(u, funwrap, 'stepInit', small);
         ixs_1 = uopt(opt.cyclical) == 1;
         ixs_0 = uopt(opt.cyclical) == 0;

      end
   end
   


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

   %dofs = ~extra.disc.isdirdofs; %% exclude dirichlet nodes

   dd = dd';
   %dd = dd(dofs);
   
   [val, oval_du, oval_dd] = obj_fun(value(u), dd(:));
   
   %% use adjoint to compute gradient
   lambda = -extra.A \ oval_dd; % A symmetric, so no transpose necessary
   
   dAdu_dd = 0; % @@ will change when including stiffness params. dependence on u
   dbdu = extra.rhs.jac{1};
   dsys_du = dAdu_dd - dbdu;
   
   grad = oval_du + dsys_du' * lambda;
   %grad = grad';
   %grad = oval_du + lambda' * dsys_du;
   
   % invert signs, since the unitBoxBFGS routine maximizes rather than
   % minimizes
   val = -val;
   grad = -grad;
end


   
   %function [v, u, history] = unitBoxBFGS(u0, f, varargin)   
%function [uu, extra] = VEM_linElast_AD(G, C, el_bc, load, varargin)
   
