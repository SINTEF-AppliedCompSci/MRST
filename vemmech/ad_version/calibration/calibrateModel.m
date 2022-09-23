function res = calibrateModel(G, props, bcond, matcond, indata, varargin)
% Calibrate boundary conditions and material parameters for the model based
% on a set of input data.
%
% SYNOPSIS:
%   res = calibrateModel(G, props, bcond, matcond, indata)
%
% PARAMETERS:
%   G        - The model grid to calibrate
%   props    - a priori material properties, being a structure with the
%              fields 'K' (bulk modulus), 'G' (shear modulus) and 'rho'
%              (material density).  Each field can be either a single value
%              (for all cells), or a vector with one entry per layer, or one
%              entry per cell.
%   bcond    - (m x 2) matrix, where 'm' is either 3 or 5.  The rows
%              represent (with the last two being optional):
%              1) orientation1 apriori, orientation1 mean_dev;
%              2) strain1 apriori,      mean_dev;
%              3) strain2 apriori,      mean_dev;
%              4) dstrain1/dz apriori,  mean_dev;  (optional)
%              5) dstrain2/dz apriori,  mean_dev;  (optional)
%              Any row with a zero mean_dev is considered an imposed
%              condition, and will not be optimized for.
%              If 'bcond' is an empty matrix, there will be no calibration on
%              strain (i.e. caused by far-field tectonic stress), which will
%              remain 0.
%   matcond  - Specify n layer, apriori multipliers and mean_dev for material
%              properties.  Structured as a matrix  with n lines and 4
%              columns, where line 'i' specifies:
%
%              [ix, K_i mean_dev, G_i mean_dev, rho_i mean_dev]
%   
%              Here, 'ix' specify the uppermost layer of cells in the grid
%              belonging to the zone covered by the multiplier
%
%              If the matrix is empty, no optimization will be done on
%              elastic parameters.
%      
%              If a mean_dev is nonpositive, no calibration will be carried
%              out for that entry.
% 
%   indata   - Indata given as a collection of 3xn matrices in a struct with
%              the following fields 'orientation', 'Smax', 'Smin' and 'Sv'.
%              Each field is a 3xn matrix where the first column indicate
%              cell index, the second indicate the actual measured value, and
%              the third indicate the associated mean deviation.
%
% RETURNS:
%   res - Structure with the following fields:
%              - bcond_opt   : vector with optimal values for boundary conditions
%              - matcond_opt : vector with optimal values for the material
%                              condition multipliers
%              - foptval     : value of the cost function at the optimum
%              - history     : optimization history
%              - uu_opt      : displacement field at optimum
%              - extra       : discretization structure (needed to compute
%                              e.g. stresses from the displacements)

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

   mrstModule add ad-core linearsolvers optimization vemmech
   
   opt.ignore_prior = false;
   opt.ignore_bnd_prior = false; % ignore prior on boundary condition
                                 % (whether or not prior is ignored for other controls)
   opt.scaling = []; % scaling objective function (should be chosen to make
                     % objective function values on the order of unity)
   opt.produce_syndata_for = []; % provide control vector, return case
   opt.gradTol = 1e-6; 
   opt.objChangeTol = 1e-7; % 1e-10(much tighter then unitBoxBFGS default) too tight?
   opt.maxIt = 100; % maximum number of directional searches in optimization
   opt.lineSearchMaxIt = 5; % maximum number of evaluations in line search
   opt.backup_file = [];
   opt.wolfe1 = 1e-3;
   opt.wolfe2 = 0.9; 
   opt.B_scale = 1; %bfgs_step = 1; % chose scaling of initial (diagonal) matrix for quasi-hessian
   opt.u_start = [];  % initial guess for start parameters
   opt.Hi = [];       % initial guess for quasi-hessian to BFGS routine
   opt.cellcount_itersolve = 50000; % how many cells in a grid before switching from
                                    % default Matlab solver to iterative
                                    % solver
   opt.model_sigma = zeros(1, 4); % modelization error (standard deviations): 
                                  % [sigma_theta, sigma_smax, sigma_smin, sigma_sv]
                                  % whereas sigma_theta is absolute and
                                  % sigma_smax/smin/sv are relative
   opt.alpha_p = []; % pore pressure, multiplied with biot's parameter.  If
                     % empty, we do not consider pore pressure

   opt = merge_options(opt, varargin{:});
   
   if isempty(matcond)
      matcond = zeros(0, 4); % still empty, but with 4 columns
   end
   if isempty(bcond)
      bcond = zeros(0, 2); % still empty, but with two columns
   end
   if any(mean(G.nodes.coords(:, 1:2)) > sqrt(eps))
      warning(['Input grid was not properly centered.  Assymetric effects ' ...
               'from depth-dependent boundary conditions possible.']);
   end
   check_nonequal_init_strains(bcond); % give user a warning if strains are
                                       % equal from the start
   
   %% setting mechanics solver to be an amgsolver
   % NB: Tolerance _has_ an influence on the precision of the derivatives.
   % Optimization routine might get problems reducing derivatives when
   % magnitude becomes comparable with 'tolerance' here.
   if G.cells.num > opt.cellcount_itersolve
      mech_solver = @(A, b) callAMGCL(A, b, 'preconditioner', 'relaxation', ...
                                      'relaxation', 'gauss_seidel', ...
                                      'solver', 'cg', ...
                                      'tolerance', 1e-7, 'maxIterations', 8000);
      
      fprintf('Choosing iterative solver.\n');
   else
      mech_solver = @(A, b) A\b; % @@   
      fprintf('Choosing standard matlab backlash solver.\n');
   end

   
   %% Determine boundary conditions
   binfo = identify_boundary_elements(G);   
   [bcfun, num_bc_ctrls] = set_bcfun(G, binfo, bcond);
   
   %% Determine material conditions
   [mfuns, num_mctrls] = ...
       set_material_fun(G, props, matcond, num_bc_ctrls);
         
   mat_facs = reshape(matcond(:, 2:end), [], 1);
   %% produce synthetic data and return if requested
   % synthetic data is useful for testing the calibration
   if ~isempty(opt.produce_syndata_for)
      res = produce_syndata(G, bcfun, mfuns, opt.produce_syndata_for, mech_solver, ...
                            opt.alpha_p);
      return
   end

   %% initialize controls
   num_controls = num_bc_ctrls + sum(num_mctrls);
   u_init = zeros(num_controls, 1);
   u_sigma = ones(num_controls, 1);
   
   if opt.ignore_bnd_prior
      u_sigma(1:num_bc_ctrls) = Inf;
   end
   
      
   %% setup objective function
   fprintf('Setting up objective function.\n');
   [obj_fun, op_0, matchfun] = ...
       setup_objective_fun(G, indata, bcfun, mfuns, num_bc_ctrls, ...
                           mat_facs, u_init, u_sigma, opt.model_sigma, ...
                           opt.scaling, opt.ignore_prior, mech_solver, ...
                           opt.alpha_p);
   
   %% run optimization
   u_start = opt.u_start;
   if isempty(u_start)
      u_start = u_init;
   end
   fprintf('Starting the optimization\n');
   [foptval, uopt, history] = optimizeSR1(u_start, obj_fun, ...
                                          'backup_file', opt.backup_file, ...
                                          'delta', 0.1, ...
                                          'epsilon', opt.gradTol, ...
                                          'funval_tol', opt.objChangeTol, ...
                                          'B_scale', opt.B_scale, ...
                                          'B_init', opt.Hi, ...
                                          'maxIt', opt.maxIt); 
   
   %% Compute initial result
   E = mfuns.efun(u_init);
   nu = mfuns.nufun(u_init);
   load = mfuns.loadfun(u_init);
   bc = bcfun(u_init);
   
   [uu_init, extra_init] = VEM_linElast_AD(G, E, nu, bc, load, 'pressure', opt.alpha_p,...
                                           'linsolve', mech_solver, 'extra', op_0);
   
   %% Compute optimal result
   E = mfuns.efun(uopt);
   nu = mfuns.nufun(uopt);
   load = mfuns.loadfun(uopt);
   bc = bcfun(uopt);
   
   [uu_opt, extra_opt] = VEM_linElast_AD(G, E, nu, bc, load, 'pressure', opt.alpha_p, ...
                                         'linsolve', mech_solver, 'extra', op_0);
   
   %% Setting up return values
   res.foptval = foptval;
   res.history = history;
   res.u_opt = uopt;
   res.uu_opt = uu_opt;
   res.extra = extra_opt;
   res.bcfun = bcfun;
   res.mfuns = mfuns;
   res.u_init = u_init;
   res.uu_init = uu_init;
   res.extra_init = extra_init;
   res.obj_fun = obj_fun;
   res.matchfun = matchfun; % to compute match
   res.indata = indata;
end   

% ----------------------------------------------------------------------------
function res = produce_syndata(G, bcfun, mfuns, u, solver, alpha_p)
   
   [uu, extra] = VEM_linElast_AD(G, mfuns.efun(u), mfuns.nufun(u), ...
                                 bcfun(u), mfuns.loadfun(u), ...
                                 'pressure', alpha_p, ...
                                 'linsolve', solver);
   res.bcfun = bcfun;
   res.mfuns = mfuns;
   res.uu = uu;
   res.extra = extra;
   
end

% ----------------------------------------------------------------------------
function [objfun, op_0, matchfun] = setup_objective_fun(G, data, bcfun, mfuns, ...
                                                        num_bc_ctrls, mat_facs, ...
                                                        u_init, u_sigma, ...
                                                        model_sigma, scaling, ...
                                                        ignore_prior, solver, ...
                                                        alpha_p)
   if isempty(data.Sv)
      data.Sv = zeros(0, 3);
   end
   if isempty(data.orientation)
      data.orientation = zeros(0, 3);
   end
   if isempty(data.Smax)
      data.Smax = zeros(0, 3);
   end
   if isempty(data.Smin)
      data.Smin = zeros(0, 3);
   end
   orientation = data.orientation;
   Smax = data.Smax;
   Smin = data.Smin;
   Sv = data.Sv;
         
   d_cell =  [orientation(:, 1); Smax(:, 1); Smin(:, 1); Sv(:,1)];
   d_obs =   [orientation(:, 2); Smax(:, 2); Smin(:, 2); Sv(:,2)];
      
   d_sigma = [orientation(:, 3); Smax(:, 3); Smin(:, 3); Sv(:,3)];
   
   d_type = [1 * ones(size(orientation, 1), 1);
             2 * ones(size(Smax, 1), 1);
             3 * ones(size(Smin, 1), 1);
             4 * ones(size(Sv, 1), 1)];

   % objective fun, which may include prior term
   ofun = @(u, xx, op) gaussian_objfun(u, xx, op, G, ...
                                       d_obs, d_sigma, d_cell, ...
                                       d_type, model_sigma, bcfun, num_bc_ctrls, mat_facs, ...
                                       alpha_p, ...
                                       mfuns.efun, mfuns.nufun, u_init, ...
                                       u_sigma, scaling, ignore_prior);
   % matchfun, which only considers fit to existing data
   matchfun = @(u, xx, op) gaussian_objfun(u, xx, op, G, ...
                                           d_obs, d_sigma, d_cell, ...
                                           d_type, model_sigma, bcfun, num_bc_ctrls, mat_facs, ...
                                           alpha_p, ...
                                           mfuns.efun, mfuns.nufun, u_init, ...
                                           u_sigma, scaling, true);
   
   % remove dependence on xx
   [uu, op_0] = VEM_linElast_AD(G, mfuns.efun(u_init), mfuns.nufun(u_init), ...
                                bcfun(u_init), mfuns.loadfun(u_init), ...
                                'pressure', alpha_p, ...
                                'no_solve', true);
   objfun = @(u) eliminate_xx_dependence(u, ofun, op_0, mfuns, bcfun, G, ...
                                         solver, alpha_p, ignore_prior);
end

% ----------------------------------------------------------------------------

function [val, grad, qHess, qHessDir] = ...
       eliminate_xx_dependence(u, ofun, extra0, mfuns, bcfun, G, solver, ...
                               alpha_p, ignore_prior)
   fprintf('Entering ''eliminate_xx_dependence\n')
   u = initVariablesADI(value(u));
   bc = bcfun(u);
   E = mfuns.efun(u);
   nu = mfuns.nufun(u);
   load = mfuns.loadfun(u);
   
   fprintf('Assembling and solving linear system inside eliminate_xx_dependence\n');
   [dd, extra] = VEM_linElast_AD(G, E, nu, bc, load, ...
                                 'linsolve', solver, ...
                                 'pressure', alpha_p, ...
                                 'extra', extra0);
   dd = dd';
   fprintf('Computing objective function value...\n');
   [val, oval_du, oval_dd] = ofun(value(u), dd(:), extra);
   fprintf('Finished computing objective function value.\n');
   %% use adjoint to compute gradient
   fprintf('Computing adjoint\n');
   tic;
   lambda = -1 * solver(extra.A, full(oval_dd));
   toc;
   fprintf('Finished computing adjoint.\n');

   dAdu_dd = extra.Ax_derivs; % @@ check this
   if isa(extra.rhs, 'ADI')
      dbdu = extra.rhs.jac{1};
      dsys_du = dAdu_dd - dbdu;
   else
      dsys_du = dAdu_dd;
   end
      
   fprintf('Computing obj.function gradient...\n');
   grad = oval_du + dsys_du' * lambda;
   fprintf('Finished computing gradient.\n');
   
   %% invert signs, since the unitBoxBFGS routine maximizes rather than
   % minimizes
   % val = -1 * val;
   % grad = -1 * grad;
   % qHess = -1 * qHess; % only non-empty if explicitly requested
   fprintf('exiting eliminate_xx_dependence\n');
end

% ----------------------------------------------------------------------------
function dth = thetadiff(theta1, theta2)
   
   dth = abs(theta1 - theta2);
   
   % ensure all dth in interval [0, pi]
   while any(dth > pi)
      ixs = dth > pi;
      dth(ixs) = dth(ixs) - pi * ones(sum(ixs), 1);
   end
   
   % if angle difference is larger than pi/2, consider complementary angle
   ixs = dth > pi/2;
   dth(ixs) = pi * ones(sum(ixs), 1) - dth(ixs);
   
end


% ----------------------------------------------------------------------------

function [val, du, dx, terms, simvals, datapts, all_simvals] = ...
       gaussian_objfun(u, xx, op, G, d_obs, d_sigma, ...
                       d_cell, d_type, model_sigma, bcfun, ...
                       num_bc_ctrls, mat_facs, alpha_p, efun, ...
                       nufun, u_prior, u_sigma, ...
                       scaling, ignore_prior)

   % compute stress tensor (AD in terms of u and xx)
   stress = calculate_stress(u, xx, G, bcfun, efun, nufun, op);
   
   % compute Smax, Smin, Sv and theta
   ixs = 0:6:numel(value(stress))-1;
   [xx_ix, yy_ix, zz_ix, xy_ix] = deal(1, 2, 3, 4);
   Sv = stress(ixs + zz_ix);
   fprintf('Computing Smax, Smin and theta');
   [Smax, Smin, theta] = ...
       Sigma2SmaxSminTheta(stress(ixs + xx_ix), ...
                           stress(ixs + yy_ix), ...
                           stress(ixs + xy_ix));
   
   % correct values for pore pressure
   if ~isempty(alpha_p)
      % the computed stress tensor is in terms of effective stress (computed
      % from displacements and elastic moduli).  The measurements are in
      % terms of total stress.  In order to compare, we need to convert out
      % effective stress values to total stress, by subtracting alpha x
      % pressure
      Smax = Smax - alpha_p;
      Smin = Smin - alpha_p;
      Sv = Sv - alpha_p;
   end
   
   theta_ixs = find(d_type == 1); theta_cells = d_cell(theta_ixs);
   Smax_ixs  = find(d_type == 2); Smax_cells  = d_cell(Smax_ixs);
   Smin_ixs  = find(d_type == 3); Smin_cells  = d_cell(Smin_ixs);
   Sv_ixs    = find(d_type == 4); Sv_cells    = d_cell(Sv_ixs);

   % identify measuring points for which we are close to isotropic stress
   threshold_theta = 5e-2; % when lower than 5% differential stress, reduce
                           % impact of angle mismatch
   
   relDs = @(cells, th) abs(Smax(cells) - Smin(cells)) ./ (th * abs(Smax(cells)));
   
   theta_relDs = relDs(theta_cells, threshold_theta);     thetaIso = theta_relDs < 1;
      
   % compute contribution of theta mismatches (pi/2 - theta since d_obs is
   % given in N°E whereas angles are internally computed with reference to x-axis)
   theta_sigma = sqrt(d_sigma(theta_ixs).^2 + model_sigma(1).^2);
   
   r_theta = thetadiff(pi/2 - d_obs(theta_ixs), theta(theta_cells)) ./ theta_sigma;
   r_theta_iso = pi/2 * ones(sum(thetaIso), 1) ./ theta_sigma(thetaIso);
   
   onevec = ones(sum(thetaIso), 1);
   scaling_fun = 0.5 * (onevec - cos(theta_relDs(thetaIso)*pi));
   
   r_theta(thetaIso) =  scaling_fun .* r_theta(thetaIso) + (onevec-scaling_fun) .* r_theta_iso;
       
   
   % compute contribution of prinicpal stress mismatches
   Smax_sigma = sqrt((d_sigma(Smax_ixs) .* d_obs(Smax_ixs)).^2 + model_sigma(2).^2);
   r_Smax = (d_obs(Smax_ixs) - Smax(Smax_cells)) ./ Smax_sigma;

   Smin_sigma = sqrt((d_sigma(Smin_ixs) .* d_obs(Smin_ixs)).^2 + model_sigma(3).^2);
   r_Smin = (d_obs(Smin_ixs) - Smin(Smin_cells)) ./ Smin_sigma;
   
   Sv_sigma = sqrt((d_sigma(Sv_ixs) .* d_obs(Sv_ixs)).^2 + model_sigma(4).^2);
   r_Sv = (d_obs(Sv_ixs) - Sv(Sv_cells)) ./ Sv_sigma;

   % adding contributions from data mismatch
   match_term = sum(r_theta.^2) + ...
                sum(r_Smin.^2) + ...
                sum(r_Smax.^2) + ...
                sum(r_Sv.^2);

   num_terms = numel(d_obs);
   
   % adding contribution from prior
   u = initVariablesADI(u);
      
   r_prior = (u - u_prior(:))./ u_sigma; % correct for boundary params
      
   % for material parameters, we will have to consider log-space
   mat_ind = false(numel(value(u)), 1);
   mat_ind(num_bc_ctrls+1:end) = true;
   if any(mat_ind)
      r_prior(mat_ind) = log(1 + r_prior(mat_ind) .* mat_facs) ./ ...
          sqrt(log(1 + mat_facs.^2));
   end
      
   % 'r_prior' does not depend on xx, but make it compatible with
   % 'match_term' in order to add up below
   r_prior.jac = [r_prior.jac, {sparse(numel(value(r_prior)), ...
                                       size(match_term.jac{2}, 2))}];

   if ~ignore_prior         
      match_term = match_term + sum(r_prior.^2);
      %num_terms = num_terms + numel(value(u));
   end
   
   % computing final objective value and derivatives  
   if isempty(scaling)
      scaling = 1;
   end
   val = match_term / (num_terms * scaling);
   du = val.jac{1}(:);
   dx = val.jac{2}(:);
   val = value(val);
   
   % return all the individual terms (as AD-variables), if requested
   if nargout > 3
      %terms = [r_theta; r_Smax; r_Smin; r_Sv; r_prior] / sqrt(num_terms * %scaling);
      terms.r_theta = r_theta;
      terms.r_Smax = r_Smax;
      terms.r_Smin = r_Smin;
      terms.r_Sv = r_Sv;
      terms.r_prior = r_prior;
   end
   if nargout > 4
      % total stress (not effective stress) whenever pore pressures are involved
      simvals.Smax = Smax(Smax_cells);
      simvals.Smin = Smin(Smin_cells);
      simvals.Sv = Sv(Sv_cells);
      simvals.theta = theta(theta_cells);
   end
   if nargout > 5
      datapts.theta = [d_cell(theta_ixs), d_obs(theta_ixs), theta_sigma];
      datapts.Smax = [d_cell(Smax_ixs), d_obs(Smax_ixs), Smax_sigma];
      datapts.Smin = [Smin_cells(:), d_obs(Smin_ixs), Smin_sigma];
      datapts.Sv = [Sv_cells(:), d_obs(Sv_ixs), Sv_sigma];
   end
   if nargout > 6
      % total stress (not effective stress) whenever pore pressures are involved
      all_simvals.theta = theta;
      all_simvals.Smax = Smax;
      all_simvals.Smin = Smin;
      all_simvals.Sv = Sv;
   end
end

% ----------------------------------------------------------------------------
function check_nonequal_init_strains(bcond)

   % ensure the two strains are not equal from the start
   if ~isempty(bcond) && ... % optimize for strains
      abs(bcond(2,1) - bcond(3,1)) / mean(abs(bcond(2:3,1))) < 1e-3 && ...
      (bcond(2,2) > 0 || bcond(3,2) > 0) % at least one active strain control
      
      if bcond(1, 2) == 0
         % we are not optimizing on orientation, so nonequal strains are not
         % a concern
         return;
      end
      
      warning(['We are optimizing on orientation, but strains are set out as ' ...
             'equal from departure.  Indeterminate angles and ensuing ' ...
               'non-convergence are likely!']);
   end
end

% % ----------------------------------------------------------------------------
% function binfo = identify_boundary_elements(G)

%    sides = {'xmin', 'xmax', 'ymin', 'ymax', 'zmin', 'zmax'};
%    side_faces = cell(6, 1);
%    side_nodes = cell(6, 1);
   
%    for i = 1:numel(sides)
%       side = sides{i};
%       tmp = pside([], G, side, 1);
%       if isempty(tmp)
%          binfo = [];
%          return; % was not able to identify boundary elements
%       end
      
%       side_faces{i} = tmp.face;
%       side_nodes{i} = ...
%           unique(G.faces.nodes(mcolon(G.faces.nodePos(side_faces{i}), ...
%                                       G.faces.nodePos(side_faces{i}+1)-1)));
%    end
%    binfo.side_faces = side_faces;
%    binfo.side_nodes = side_nodes;
% end
