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

   mrstModule add optimization linearsolvers;
   
   opt.gradTol = 1e-3; %1e-3;
   opt.objChangeTol = 1e-5;%1e-5; %@@ might be too tight (much tighter than default
                            %in unitBoxBFGS)
   opt.cyclical = []; % indices of cyclical control variables
   opt.extra = []; % if discretization is precomputed, it can be passed in
                   % here to save time in VEM_linelast_AD
   [opt, ~] = merge_options(opt, varargin{:});

   if ~strcmpi(G.type{end}, 'createAugmentedGrid')
      G = createAugmentedGrid(computeGeometry(G));
   end
   
   funwrap = @(u) fun_wrapper(u, G, bcfun, cfun, loadfun, obj_fun, opt.extra);
                              
   
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

      amgsolver = @(A, b) callAMGCL(A, b, 'relaxation', 'chebyshev', 'solver', 'cg', ...
                                    'tolerance', 2e-6, 'maxIterations', 2000);
      if nargout == 4
         uu_opt = VEM_linElast_AD(G, C, bc, load, 'linsolve', amgsolver, ...
                                  'extra', opt.extra);
      else
         [uu_opt, extra] = VEM_linElast_AD(G, C, bc, load, 'linsolve', amgsolver, ...
                                           'extra', opt.extra);
      end
   end
end

function [val, grad] = fun_wrapper(u, G, bcfun, cfun, loadfun, obj_fun, extra)

   fprintf('Calling fun_wrapper\n');
   u = initVariablesADI(u);
   
   bc = bcfun(u);
   C = cfun(u);
   load = loadfun(u);
   
   amgsolver = @(A, b) callAMGCL(A, b, 'relaxation', 'chebyshev', 'solver', 'cg', ...
                              'tolerance', 2e-6, 'maxIterations', 2000);
   
   [dd, extra] = VEM_linElast_AD(G, C, bc, load, ...
                                 'linsolve', amgsolver, 'extra', extra);

   %dofs = ~extra.disc.isdirdofs; %% exclude dirichlet nodes

   dd = dd';
   %dd = dd(dofs);
   
   [val, oval_du, oval_dd] = obj_fun(value(u), dd(:));
   
   %% use adjoint to compute gradient
   %   keyboard;
   tic; fprintf('Solving for adjoint.\n');
   lambda = amgsolver(extra.A, -full(oval_dd));
   %lambda = -extra.A \ oval_dd; % A symmetric, so no transpose necessary
   toc;
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