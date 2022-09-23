function [foptval, uopt, history, uu_opt, extra] = ...
       optimize_mech(u, G, bcfun, efun, nufun, loadfun, obj_fun, varargin)
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
%   efun    - function that generates Young's modulus from the control
%             parameters, for all cells in the grid
%   nufun   - function that generates Poisson's parameter from the control
%             parameters, for all cells in the grid
%   loadfun - function that generates the cellwise load term (AD) from the
%             control parameters
%   obj_fun - objective function to minimize.  Should take the vector of
%             control variables as its first argument, the computed
%             displacements as its second argument, and the VEM discretization
%             ('extra') as its third argument (to avoid need of recomputing
%             it), and return the following:
%             + val - the objective value
%             + du  - partial derivatives of the objective value wrt. control
%                     variables
%             + dd  - partial derivatives of the objective value wrt. 
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

%{
Copyright 2009-2022 SINTEF Digital, Mathematics & Cybernetics.

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

   mrstModule add optimization linearsolvers;
   
   opt.gradTol = 1e-9; %1e-3;
   opt.objChangeTol = 1e-10;%1e-10;%1e-5; %@@ might be too tight (much tighter than default
                            %in unitBoxBFGS)
   opt.extra = []; % if discretization is precomputed, it can be passed in
                   % here to save time in VEM_linelast_AD
   opt.background_forces = []; % background forces acting on selected nodes,
                               % and which are not adjusted during optimization
   [opt, ~] = merge_options(opt, varargin{:});

   if ~strcmpi(G.type{end}, 'createAugmentedGrid')
      G = createAugmentedGrid(computeGeometry(G));
   end
   
   funwrap = @(u) fun_wrapper(u, G, bcfun, efun, nufun, loadfun, obj_fun, opt.extra, ...
                              opt.background_forces);

   
   [foptval, uopt, history] = optimizeSR1(u, funwrap, ...
                                          'epsilon', opt.gradTol, ...
                                          'funval_tol', opt.objChangeTol);
                                          

   %% compute additional information if requested
   if nargout > 3
      E = efun(uopt);
      nu = nufun(uopt);
      bc = bcfun(uopt);
      load = loadfun(uopt);

      amgsolver = @(A, b) callAMGCL(A, b, 'relaxation', 'chebyshev', 'solver', 'cg', ...
                                    'tolerance', 2e-6, 'maxIterations', 8000);
      if nargout == 4
         uu_opt = VEM_linElast_AD(G, E, nu, bc, load, 'linsolve', amgsolver, ...
                                  'extra', opt.extra, 'background_forces', opt.background_forces);
      else
         [uu_opt, extra] = VEM_linElast_AD(G, E, nu, bc, load, 'linsolve', amgsolver, ...
                                           'extra', opt.extra, 'background_forces', ...
                                           opt.background_forces);
      end
   end
end


% ----------------------------------------------------------------------------
function [val, grad] = fun_wrapper(u, G, bcfun, efun, nufun, loadfun, ...
                                   obj_fun, extra, background_forces)

   fprintf('Calling fun_wrapper\n');
   u = initVariablesADI(value(u));
   
   bc = bcfun(u);
   E = efun(u);
   nu = nufun(u);
   load = loadfun(u);

   amgsolver = @(A, b) callAMGCL(A, b, 'relaxation', 'chebyshev', 'solver', 'cg', ...
                              'tolerance', 2e-6, 'maxIterations', 8000);
   
   [dd, extra] = VEM_linElast_AD(G, E, nu, bc, load, ...
                                 'linsolve', amgsolver, ...
                                 'extra', extra, ...
                                 'background_forces', background_forces);
                                 %'force_method', 'cell_force', ...
   dd = dd';

   fprintf('Calling obj_fun\n');
   [val, oval_du, oval_dd] = obj_fun(value(u), dd(:), extra);

   
   %% use adjoint to compute gradient
   tic; fprintf('Solving for adjoint.\n');
   lambda = -1 * amgsolver(extra.A, full(oval_dd));
   %lambda = -extra.A \ oval_dd; % A symmetric, so no transpose necessary
   toc;
   dAdu_dd = extra.Ax_derivs; % @@ check this
   if isa(extra.rhs, 'ADI')
      dbdu = extra.rhs.jac{1};
      dsys_du = dAdu_dd - dbdu;
   else
      dsys_du = dAdu_dd;
   end
      
   grad = oval_du + dsys_du' * lambda;

   % % take care of possible scaling issue if derivatives become unreasonable
   % % large (may typically happen at very first step if initial guess is
   % 'bad')
   
   if norm(grad) > 1./sqrt(eps)
      grad = grad / sqrt(norm(grad));
   end
   
   
end
