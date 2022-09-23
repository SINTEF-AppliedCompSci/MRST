function [objfun, op_0, matchfun, ofun] = setup_objective_fun(G, data, bcfun, mfuns, ...
                                                                 num_bc_ctrls, mat_facs, ...
                                                                 u_init, u_sigma, ...
                                                                 model_sigma, scaling, ...
                                                                 ignore_prior, solver, ...
                                                                 alpha_p)
%Undocumented Utility Function

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
      
   % compute contribution of theta mismatches 
   theta_sigma = sqrt(d_sigma(theta_ixs).^2 + model_sigma(1).^2);
   
   r_theta = thetadiff(d_obs(theta_ixs), theta(theta_cells)) ./ theta_sigma;
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
