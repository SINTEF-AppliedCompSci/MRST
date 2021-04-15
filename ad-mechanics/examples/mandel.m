%% Purpose of this script
% This script is the accompanying software for the third example discussed in
% Chapter 14, "A Brief Introduction to Poroelasticity and Simulation of Coupled
% Geomechanics and Flow in MRST", of the book "Advanced Modelling with the
% MATLAB Reservoir Simulation Toolbox (MRST)". The script sets up and run a
% simulation of Mandel's problem, first described by Jean Mandel in:
% 
% Mandel, J. (1953) "Consolidation des sols (étude mathématique)"
% Géotechnique 3, 287-299.
%

%% Include necessary modules
mrstModule add ad-mechanics ad-core ad-props ad-blackoil vemmech
gravity off; % we do not consider gravity in this example

%% define grid
L = 10 * meter;
H = 1 * meter;
physdim = [L, H];
resolution = [50, 10];

G = cartGrid(resolution, physdim);
G = computeGeometry(G);

%% Plot grid
figure; 
plotGrid(G);
set(gcf, 'position', [1, 1, 900, 400]);
set(gca, 'fontsize', 15);

%% Setup rock parameters

% flow parameters
perm = 200 * milli * darcy;
poro = 0.1;

% elastic parameters
E = 0.1 * giga * Pascal;   % Young's modulus
nu = 0.2;                  % Poisson's parameter
alpha = 1;                 % Biot Willis coefficient

rock.poro  = poro  * ones(G.cells.num, 1);
rock.perm  = perm  * ones(G.cells.num, 1);
rock.alpha = alpha * ones(G.cells.num, 1);

top_press = 1 * mega * Pascal;

% boundary conditions will be set in the 'mech_problem' object further down
mech_problem.E = E * ones(G.cells.num, 1);
mech_problem.nu = nu * ones(G.cells.num, 1);
mech_problem.load = @(x) 0*x;  % no body force applied

%% setup fluid parameters

muW = 0.89 * milli * Pascal / second; % fluid viscosity
rhoW = 1000 * kilogram / meter^3;     % fluid density
cW = 1.0e-10 * Pascal^-1;             % fluid compressibility
pRef = 0;                             % reference pressure (zero)

% Conceptually, the fluid object created below treats fluid compressibility as
% an exponential function, but with physically realistic values for 'cW', we
% remain within the regime of linear compressibility to an excellent approximation.
fluid = initSimpleADIFluid('phases' , 'W'  , ... % only water phase
                           'mu'     , muW  , ...
                           'rho'    , rhoW , ...
                           'c'      , cW   , ...
                           'pRef'   , pRef);

%% Set boundary conditions

% define lambda function to identify nodes for a set of faces
facenodes = ...
    @(f) unique(G.faces.nodes(mcolon(G.faces.nodePos(f), ...
                                     G.faces.nodePos(f+1)-1)));

% identify bottom faces and nodes
bfaces = find(G.faces.centroids(:,2) == min(G.faces.centroids(:,2)));
bnodes = facenodes(bfaces); % bottom nodes
nbn = numel(bnodes); % number of bottom nodes (and top nodes)

% identify top faces and nodes
tfaces = find(G.faces.centroids(:,2) == max(G.faces.centroids(:,2)));
tnodes = facenodes(tfaces);

% identify left boundary faces and nodes
xmin_faces = find(G.faces.centroids(:,1) == 0);
lnodes = facenodes(xmin_faces);
nln = numel(lnodes); % number of nodes on left side

% identify right boundary faces (no need for nodes)
xmax_faces = find(G.faces.centroids(:,1) == max(G.faces.centroids(:, 1)));

% define flow boundary conditions (fixed pressure on right boundary, no-flow
% otherwise)
bc = addBC([], xmax_faces, 'pressure', pRef);

% setup displacement boundary conditions (actual displacement value for top
% nodes to be determined later)
disp_bc = struct('nodes', [1:G.nodes.num]', ...
                 'uu',  zeros(size(G.nodes.coords)), ...
                 'mask', false(G.nodes.num, 2));
disp_bc.mask(tnodes, 2) = true;
disp_bc.mask(bnodes, 2) = true;
disp_bc.mask(lnodes, 1) = true;

%% setup initial state and determine time scale for simulation

% initial state
mech_problem.el_bc = struct('disp_bc', disp_bc, 'force_bc', []);
model = MechWaterModel(G, rock, fluid, mech_problem);
mech_unknowns = ~model.mechModel.operators.isdirdofs;

initState.pressure = pRef * ones(G.cells.num, 1);
initState.xd = zeros(nnz(mech_unknowns), 1);
initState = addDerivedQuantities(model.mechModel, initState);

%%compute characteristic time scale (L^2/c)

% we need to compute S (uniaxial specific storage).  One quick way to 
% do this is to call 'poroParams' with the parameters we already have.
params = poroParams(poro, true, 'E', E, 'alpha', alpha, 'nu', nu, 'K_f', 1/cW);

c = perm / (muW * params.S);
Tchar = L^2 / c; % charateristic time

%% Run simulation

% identify top cell indices
tcells = sum(G.faces.neighbors(tfaces, :), 2);

tsteps = 200;
duration = 2 * Tchar * second; % total duration (two times characteristic time)

dy_init = -1e-3; % initial guess for vertical displacement
state0 = initState;
states = {};

schedule.step.val = duration/tsteps; % a single step
schedule.step.control = 1;
schedule.control = struct('W', [], 'bc', bc);

for step = 1:tsteps 
   % compute state for next timestep
   fprintf('\nSolving timestep %i\n\n', step);
   while (true)
      
      disp_bc.uu(tnodes, 2) = dy_init;
      
      % combine displacement and force boundary conditions
      mech_problem.el_bc = struct('disp_bc', disp_bc, 'force_bc', []);
      
      % model
      model = MechWaterModel(G, rock, fluid, mech_problem);

      % define an initial guess, to ensure the correct values for imposed
      % displacements are respected when the timestep is simulated below
      initGuessState = addDerivedQuantities(model.mechModel, state0);
      
      % simulate a single timestep
      [~, state] = simulateScheduleAD(state0, model, schedule, ...
                                      'initialGuess', initGuessState);

      % comparing simulated mean stress with the target      
      tstress_eff = state{1}.stress(tcells, 2);
      tpress  = state{1}.pressure(tcells);
      
      tforce_sim = mean(-tstress_eff + alpha * tpress);

      if abs((top_press - tforce_sim)/top_press) < 1e-4
         % the simulated pressure was close enough to the target.  Save the
         % current simulated timestep, and proceed to next step.
         break;
      else
         % adjust displacement boundary condition and try again
         dy_init = dy_init * top_press / tforce_sim;
      end
   end
   states = [states, state];
   state0 = state{1};
end


%% Visualize intermediary result
figure;
plotCellData(G, states{100}.pressure); colorbar
set(gcf, 'position', [1, 1, 900, 400]);
set(gca, 'fontsize', 15);

%% Plot pressure evolution over time
figure; hold on;
dimless_time = [0.01, 0.1, 0.5, 1.0, 2.0];
styles = {'b-', 'rx-', 'mo-', 'k--', 'gs'};
abscissa = [G.cells.centroids(tcells, 1); L] / L;
ix = 1;
for tau = dimless_time
   tstep_ix = ceil(tau * tsteps/2);  
   profile = [states{tstep_ix}.pressure(tcells); 0] / (top_press / 2);
   plot(abscissa, profile, styles{ix}, 'linewidth', 1.5);
   ix = ix+1;
end
legend('\tau=0.01', '\tau=0.1', '\tau=0.5', '\tau=1.0', '\tau=2.0');
xlabel('x/L');
ylabel('2p/p_0');
set(gca, 'fontsize', 15);
set(gcf, 'position', [0, 0, 850, 500]);
axis tight;

% nb: above, top_stress is divided by 2 to give 2D mean stress.  Skempton's
% coefficient (here 1) gives the ratio of p to mean stress (not sigma_zz)

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
