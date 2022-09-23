mrstModule add vemmech ad-core linearsolvers mrst-gui optimization
gravity on;

%% Construct test grid

gridres = [10, 10, 10];

griddims = [3, 3, 1] * kilo * meter;
props.young = 5e9;
props.rho = 1500; 
props.poisson = 0.25;

G = cartGrid(gridres, griddims);
G = createAugmentedGrid(computeGeometry(G));

% expand properties to every grid cell
for f = {'young', 'rho', 'poisson'}
    props.(f{:}) = props.(f{:}) * ones(G.cells.num, 1);
end

% compute elasticity tensor
props.C = Enu2C(props.young, props.poisson, G);


%% Specify control variables and boundary condition function

binfo = identify_boundary_elements(G);

theta_cond = [0, pi/4]; % mean and rate of change with control param.
strain_cond1 = [0, 1e-4]; % mean and rate of change
strain_cond2 = [-2e-4, 1e-4]; % as above, but slightly different to avoid
                               % initial isotropy
% strain_cond1 = [-2e-5, 1e-5]; % mean and rate of change
% strain_cond2 = [-1e-5, 1e-5]; % as above, but slightly different to avoid
%                                % initial isotropy
bcond = vertcat(theta_cond, strain_cond1, strain_cond2);
num_bc_ctrls = size(bcond, 1);

% construct the boundary condition function, which depend on 3 control
% variables (theta, strain1 and strain2)
[bcfun, num_bc_ctrls] = set_bcfun(G, binfo, bcond);

num_layers = 3;
layers_start_z_ix = [1; 4; 7]; % update this if grid z-resolution changes!

% layers_start_z_ix = 1:floor(gridres(3)/num_layers):gridres(3);
% layers_start_z_ix = layers_start_z_ix(1:end-1);

change_rate_material = 1; % how parameter property values (bulk mod.,
                          % shear.mod and density)
                          % change (relative) with control parameters

matcond = [layers_start_z_ix, ones(num_layers, 3) * change_rate_material];

% construct the function specifying layer-wise material properties
[mfuns, num_mctrls] = set_material_fun(G, props, matcond, num_bc_ctrls);

%% Define the synthetic measurements for matching

% "unknown" control vector we are seeking to match
u_target = [1/6, 0.4, 0.6, ... % bc
            0.2, 0.5, 0.3, ... % E/K
            0.7, 0.1, 0.2, ... % nu/G
            0.2, 0.5, 0.4]'; % load

[uu_target, op] = VEM_linElast_AD(G, mfuns.efun(u_target), ...
                                  mfuns.nufun(u_target), ...
                                  bcfun(u_target), ...
                                  mfuns.loadfun(u_target));
                              
stress = calculateStressVEM(G, uu_target, op);
[Smax, Smin, theta] = Sigma2SmaxSminTheta(stress(:,1), ...
                                          stress(:,2), ...
                                          stress(:,4));
Sv = stress(:,3);

ncells = G.cells.num;
weight = 1; % we use same weight for all data points

% we use all computed values as input to the optimization.  In principle, we
% could choose a smaller subset, and different subsets for the different
% datatypes
subset = (1:ncells)';

data.Sv          = [subset, Sv(subset),    weight * ones(size(subset))];
data.orientation = [subset, theta(subset), weight * ones(size(subset))];
data.Smax        = [subset, Smax(subset),  weight * ones(size(subset))];
data.Smin        = [subset, Smin(subset),  weight * ones(size(subset))];

%% Run optimization

% setup objective function
mat_facs = reshape(matcond(:, 2:end), [], 1);
u_init   = zeros(size(u_target)); % all zeros;
u_sigma  = Inf(size(u_target)); % we do not consider uncertainty
model_sigma = zeros(1, 4); % we do not consider uncertainty in model either
scaling = [];

u_init(2) = 0.4;
u_sigma(2) = 1;
ignore_prior = true; % we pay no attention to the prior when determining
                     % optimum
solver = @(A, b) A\b;

obj_fun = setup_objective_fun(G, data, bcfun, mfuns, num_bc_ctrls, mat_facs, ...
                              u_init, u_sigma, model_sigma, scaling, ignore_prior, ...
                              solver, []);
                                 

% running optimization
[foptval, uopt, history] = optimizeSR1(u_init, obj_fun, 'delta', 0.1, ...
                                       'B_scale', 0.01, ...
                                       'r', 1e-3, ...
                                       'rat_lim', 0.75, ... 
                                       'delta_fac', 0.8 , ... 
                                       'epsilon', 1e-16, ... % grad. tolerance
                                       'funval_tol', 1e-16); % obj. change tol

%% Testing result

[uu_opt, op] = VEM_linElast_AD(G, mfuns.efun(uopt), ...
                                  mfuns.nufun(uopt), ...
                                  bcfun(uopt), ...
                                  mfuns.loadfun(uopt));
                              
stress_opt = calculateStressVEM(G, uu_opt, op);
[Smax_opt, Smin_opt, theta_opt] = Sigma2SmaxSminTheta(stress_opt(:,1), ...
                                                      stress_opt(:,2), ...
                                                      stress_opt(:,4));
Sv_opt = stress_opt(:,3);

% Comparing with target
optvals = {theta_opt, Smax_opt, Smin_opt, Sv_opt};
targetvals = {theta, Smax, Smin, Sv};
labels = {'theta', 'Smax', 'Smin', 'Sv'};

for i=1:4
    fprintf(['Relative difference in mean, ', labels{i}, ': %i\n'], ...
            abs(mean(targetvals{i}) - mean(optvals{i}))/mean(targetvals{i}));

    fprintf(['Max. relative difference, ', labels{i}, ': %i\n'], ...
            max(abs((targetvals{i}-optvals{i})./targetvals{i})));
end

% Comparing control vectors
fprintf(['Boundary conditions (left column: target, right column: search ' ...
         'result):\n']);
[u_target(1:3), uopt(1:3)]
fprintf(['Material properties:\n']);
[u_target(4:end), uopt(4:end)]

fprintf(['We can see that even though the stress field closely corresponds to ' ...
         'the target, the controls do not always.  The problem is underdetermined ' ...
         'without a prior.  We will re-run the optimization while locking one ' ...
         'of the boundary displacements to the target value.\n']);

%% Re-run optimization with prior on a single boundary variable

u_init(2) = 0.4;
u_sigma(2) = 1e-2;
ignore_prior = false; % we pay no attention to the prior when determining
                     % optimum
solver = @(A, b) A\b;

obj_fun = setup_objective_fun(G, data, bcfun, mfuns, num_bc_ctrls, mat_facs, ...
                              u_init, u_sigma, model_sigma, scaling, ignore_prior, ...
                              solver, []);
                                 

% running optimization
[foptval, uopt, history] = optimizeSR1(u_init, obj_fun, 'delta', 1, ...
                                       'B_scale', 1, ...
                                       'r', 1e-3, ...
                                       'rat_lim', 0.75, ... 
                                       'delta_fac', 0.8 , ... 
                                       'epsilon', 1e-16, ... % grad. tolerance
                                       'funval_tol', 1e-16); % obj. change tol

[uu_opt, op] = VEM_linElast_AD(G, mfuns.efun(uopt), ...
                                  mfuns.nufun(uopt), ...
                                  bcfun(uopt), ...
                                  mfuns.loadfun(uopt));
                              
stress_opt = calculateStressVEM(G, uu_opt, op);
[Smax_opt, Smin_opt, theta_opt] = Sigma2SmaxSminTheta(stress_opt(:,1), ...
                                                      stress_opt(:,2), ...
                                                      stress_opt(:,4));
Sv_opt = stress_opt(:,3);

% Comparing with target
optvals = {theta_opt, Smax_opt, Smin_opt, Sv_opt};
targetvals = {theta, Smax, Smin, Sv};
labels = {'theta', 'Smax', 'Smin', 'Sv'};

for i=1:4
    fprintf(['Relative difference in mean, ', labels{i}, ': %i\n'], ...
            abs(mean(targetvals{i}) - mean(optvals{i}))/mean(targetvals{i}));

    fprintf(['Max. relative difference, ', labels{i}, ': %i\n'], ...
            max(abs((targetvals{i}-optvals{i})./targetvals{i})));
end

% Comparing control vectors
fprintf(['Boundary conditions (left column: target, right column: search ' ...
         'result):\n']);
[u_target(1:3), uopt(1:3)]
fprintf(['Material properties:\n']);
[u_target(4:end), uopt(4:end)]

fprintf('We note that the correct target controls are now identified.\n');

%% Re-run, to demonstrate use of optimize_mech function

% It is also possible to use the optimize_mech function for optimization,
% which may be a bit easier if you want to implement your own objective
% function and do not want to handle the computation of adjoint-based
% gradients yourself.  Note that in the call to setup_objective_function
% below, we only return the objective function as a function of control
% variables _and_ displacements here.

[~, ~, ~, ofun] = setup_objective_fun(G, data, bcfun, mfuns, num_bc_ctrls, mat_facs, ...
                                      u_init, u_sigma, model_sigma, scaling, ignore_prior, ...
                                      solver, []);
                                 
[foptval, uopt, history] = ...
    optimize_mech(u_init, G, bcfun, mfuns.efun, mfuns.nufun, mfuns.loadfun, ofun);

%% Copyright Notice
%
% <html>
% <p><font size="-1">
% Copyright 2009-2022 SINTEF Digital, Mathematics & Cybernetics.
% </font></p>
% <p><font size="-1">
% This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).
% </font></p>
% <p><font size="-1">
% MRST is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% </font></p>
% <p><font size="-1">
% MRST is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% </font></p>
% <p><font size="-1">
% You should have received a copy of the GNU General Public License
% along with MRST.  If not, see
% <a href="http://www.gnu.org/licenses/">http://www.gnu.org/licenses</a>.
% </font></p>
% </html>
