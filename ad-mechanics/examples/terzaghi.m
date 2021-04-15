%% Purpose of this script
% This script is the accompanying software for the first and second examples
% discussed in Chapter 14, "A Brief Introduction to Poroelasticity and
% Simulation of Coupled Geomechanics and Flow in MRST", of the book "Advanced
% Modelling with the MATLAB Reservoir Simulation Toolbox (MRST)". The script
% sets up and run a simulation of Terzaghi's problem, first described by 
% Karl Terzaghi in:
%
% Terzaghi K. (1925) Erdbaumechanik auf Bodenphysikalischer Grundlage, Deuticke:
% Leipzig

%% Load necessary modules
mrstModule add vemmech
gravity off; % set to gravity 'on' to include gravity.  Be aware that the
             % analytic solutions of the model problems below assume no gravity.

%% First case: compression of a dry sample
% This is a purely linear elastic problem (no poroelasticity)

% Define vertical cylinder
P = [];
layers = 20;
for r = linspace(0.2, 1, 5)
   [x, y, ~] = cylinder(r^1.5, 16);
   P = [P [x(1,:); y(1,:)]];
end
P = unique([P'; 0, 0], 'rows');
aG = pebi(triangleGrid(P));
G = makeLayeredGrid(aG, layers);
G = computeGeometry(G); % compute basic geometric information
G = createAugmentedGrid(G); % compute additional geometric information needed
                            % by the VEM mechanics code

% plot the cylinder
figure
plotGrid(G, 'FaceColor', [.8, .8, .8]); view(33, 10); axis tight
set(gcf, 'color', 'white', 'position', [0 0 500 800])
set(gca, 'fontsize', 15)


% define basic parameters
density = 2000 * kilogram / (meter^3); % unless gravity is involved, this
                                       % value will not matter
Nc = G.cells.num;
E = 5 * giga * Pascal; % Young's modulus
nu = 0.3;  % Poisson's parameter

% identify bottom notes, side nodes, top nodes and top faces
bottom_nodes = find(G.nodes.coords(:,3) == max(G.nodes.coords(:,3)));
top_nodes    = find(G.nodes.coords(:,3) == 0);
top_faces    = find(G.faces.centroids(:,3) == 0);

% identify innermost bottom nodes
bottom_innermost = ...
    find(sqrt(sum(G.nodes.coords(bottom_nodes, 1:2).^2, 2)) < 0.1);

% zero vertical displacement for bottom nodes (and zero displacement for
% innermost node, to anchor the problem)
el_bc.disp_bc.nodes = bottom_nodes;
el_bc.disp_bc.uu = repmat([0, 0, 0], numel(bottom_nodes), 1);
el_bc.disp_bc.mask = repmat([false, false, true], numel(bottom_nodes), 1);
el_bc.disp_bc.mask(bottom_innermost, 1:2) = true;

% force applied at top boundary
top_force = 1e7 * Pascal;
el_bc.force_bc.faces = top_faces;
el_bc.force_bc.force = repmat([0, 0, top_force], numel(top_faces), 1);

% in the default case, we ignore gravity, but we need to specify a load
% function nevertheless.
load = @(x) repmat(density * gravity(), size(x, 1), 1);

% solving the linear elastic system 
C = Enu2C(E * ones(Nc, 1), nu * ones(Nc, 1), G);
uu = VEM_linElast(G, C, el_bc, load);

% plot result
figure
plotNodeDataDeformed(G, sqrt(sum(uu.^2, 2)), uu * 100); 
view(0, 0); colorbar; axis tight;
set(gcf, 'color', 'white', 'position', [0 0 500 800])
set(gca, 'fontsize', 15)
set(gca, 'zlim', [0, 20])


% Verify poisson ratio

% axial displacement per unit:
L = max(G.nodes.coords(:,3)) - min(G.nodes.coords(:,3)); % length
R = max(G.nodes.coords(:,1)); % should equal 1

axial_strain = uu(top_nodes(1), 3) / L;
radial_strain = max(uu(:,1)) / R;

measured_nu = radial_strain / axial_strain;
relative_err = (measured_nu - nu)/nu;
fprintf('Real nu: %1.5f.  Measured nu: %1.5f\n', nu, measured_nu);
fprintf('Relative error: %1.2e\n', relative_err);

% re-running the experiment, imposing roller boundary conditions on lateral
% boundary. 
side_nodes = find(sqrt(sum(G.nodes.coords(:,1:2).^2, 2)) > 0.9);
side_nodes = setdiff(side_nodes, bottom_nodes);
[Nb, Ns] = deal(numel(bottom_nodes), numel(side_nodes));

el_bc.disp_bc.nodes = [bottom_nodes; side_nodes];
el_bc.disp_bc.uu = repmat([0, 0, 0], Nb + Ns, 1);
el_bc.disp_bc.mask = [repmat([true, true, true], Nb, 1); ... % locked btm.
                      repmat([true, true, false], Ns, 1)]; % roller side bnd.
                      
% solving the linear elastic system
uu = VEM_linElast(G, C, el_bc, load);

% plot result
figure
plotNodeDataDeformed(G, sqrt(sum(uu.^2, 2)), uu * 100); 
view(0, 0); colorbar; axis tight;
set(gcf, 'color', 'white', 'position', [0 0 500 800])
set(gca, 'fontsize', 15)
set(gca, 'zlim', [0, 20])

K = E / (3 * (1-2*nu)); % compute bulk modulus from E and nu
Kv = 3 * K * (1-nu)/(1+nu); % compute vertical incompressibility

axial_strain = uu(top_nodes(1), 3) / L;
predicted_strain = top_force / Kv;
relative_err = (axial_strain - predicted_strain)/predicted_strain;
fprintf('Predicted strain: %1.5f.  Measured strain: %1.5f\n', ...
        predicted_strain, axial_strain);
fprintf('Relative error: %1.2e\n', relative_err);


%% Terzhagi's problem
mrstModule add ad-mechanics ad-core ad-props ad-blackoil

% model
perm = 300 * milli * darcy; % permeability
poro = 1/4;                 % porosity 
pRef = 0;                   % ref. pressure for fluid density
alpha = 0.9;                % Biot Willis coefficient
Kf = 1.96 * giga * Pascal;  % fluid incompressibility

% defining poroelastic model (here using the fully coupled MechWaterModel)

rock = struct('perm', perm * ones(Nc, 1), ...
              'poro', poro * ones(Nc, 1), ...
              'alpha', alpha * ones(Nc, 1));

pvMult = (1-alpha) * (alpha-poro) / poro / K; % pore volume multiplier,
                                              % consistent with our choice of
                                              % poroelastic parameters
fluid = initSimpleADIFluid('phases', 'W', ...
                           'mu'    , 1 * centi * poise, ...
                           'rho'   , 1000 * kilogram / meter^3, ...
                           'c'     , 1 / Kf, ...
                           'cR'    , pvMult, ...
                           'pRef'  , pRef);

mech_problem = struct('E'    , E * ones(Nc, 1) , ...
                      'nu'   , nu * ones(Nc, 1), ...
                      'el_bc', el_bc, ...
                      'load' , load);

model = MechFluidFixedStressSplitModel(G, rock, fluid, mech_problem);

% initial state
num_mech_unknowns = sum(~model.mechModel.operators.isdirdofs);

initState = struct('pressure', pRef * ones(Nc, 1), ...
                   'xd', zeros(num_mech_unknowns, 1));

initState = addDerivedQuantities(model.mechModel, initState);

% schedule
bc.face  = top_faces;
bc.type  = repmat({'pressure'}, 1, numel(top_faces));
bc.value = repmat(pRef, numel(top_faces), 1); %#ok
bc.sat   = ones(numel(top_faces), 1);

% computing characteristic time

% The sequence of poroelastic parameters to compute to derive S (which we
% need to compute characteristic time) would be:

% Ks = K / (1 - alpha);  % grain compressibility
% H = K / alpha;
% S_sigma = (1 / K - 1 / Ks) + poro * (1 / Kf - 1 / Ks); % specific storage 
%                                                        % coef. at constant
%                                                        % stress
% R = 1 / S_sigma;
% B = R / H;                                    % Skempton's coefficient
% eta = (1 - 2 * nu) / (2 * (1 - nu)) * alpha;  % poroelastic stress parameter
% S = S_sigma * (1 - 4 * eta * B / 3);          % uniaxial specific storage

% To avoid going through all the steps listed above, we use the utility
% script 'poroParams', which lets us derive all the poroelastic parameters
% that can be specified based on the information available:
params = poroParams(poro, true, 'E', E, 'nu', nu, 'alpha', alpha, 'K_f', Kf);

c = perm / (fluid.muW(pRef) * params.S); % uniaxial hydraulic diffusivity
tau = c / L^2; % dimensionless time

tinystep = 1e-5;
num_tsteps = 50;
tau_steps = linspace(tinystep, 1, num_tsteps+1);
tau_steps = [0, tau_steps];
tsteps = tau_steps * L^2 / c;

schedule = ...
    struct('step', struct('val', diff(tsteps), ...
                          'control', ones(num_tsteps+1, 1)), ...
           'control', struct('W', [], 'bc', bc));

% run simulation
[~, states, report] = simulateScheduleAD(initState, model, schedule);

% compute max pressure (undrained response pressure)

% without using poroParams(), the computation of gamma would be:

% GG = E / (2 + 2 * nu);  % shear modulus
% gamma = eta / (GG * S); % loading efficiency

pmax = top_force * params.gamma;

% plot development of column pressure profile over time
z_cells = [1:G.cells.num/layers:G.cells.num]'; %#ok
z_depths = G.cells.centroids(z_cells, 3);
ixs = [0, 1, 2, 5, 20, 50] + 1;
figure; hold on
for ix = ixs
   plot([0; z_depths]/L, ...
        [0; states{ix}.pressure(z_cells)/pmax], 'linewidth', 1.5);
end
xlabel('z/L')
ylabel('\gamma \cdot \sigma_o')
legend('\tau = 10^{-5}', '\tau = 0.02', '\tau = 0.04', '\tau = 0.1', ...
       '\tau = 0.4', '\tau = 1')
set(gca, 'fontsize', 15, 'ylim', [0, 1.05])
set(gcf, 'position', [0, 0, 840, 640])

% plot vertical displacement of top surface over time
figure
w_top = [0; cellfun(@(s) s.uu(top_nodes(1), 3), states)];
plot(tsteps * c / L^2, w_top * 100, 'linewidth', 1.5);
xlabel('\tau')
ylabel('cm');
set(gca, 'fontsize', 15)
set(gcf, 'position', [0, 0, 840, 640])

% initial displacement
w_0 = w_top(2);

% theoretical initial displacement

% if we didn't make use of poroParams(), the computation would be:

% S_epsilon = S_sigma - alpha^2/K;
% Kvu = alpha / (gamma * S_epsilon);
% w_0_theoretical = top_force * L / Kvu;

w_0_theoretical = top_force * L / params.K_vu;

rel_err = (w_0 - w_0_theoretical)/w_0;

fprintf('Expected init. displacement: %1.5f m.\n', w_0_theoretical)
fprintf('Measured displacement: %1.5f m\n', w_0);
fprintf('Relative error: %1.2e\n', rel_err);

%% Copyright Notice
%
% <html>
% <p><font size="-1">
% Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.
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
