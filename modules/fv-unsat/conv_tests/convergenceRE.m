function [errorPsi, errorFlux] = convergenceRE(cells, timeLevels)
% Convergence of Richards' equation (Numerical Example 1 from the Chapter)
%
% SYNOPSIS:
%   function [errorPsi, errorFlux] = convergenceRE(cells, timeLevels)
%
% PARAMETERS:
%   cells       - Scalar, number of cells = nx = ny
%   timeLevels  - Scalar, number of discrete time levels
%
% RETURNS:
%   errorPsi    - Scalar, L2-like discrete error for the pressure
%   errorFlux   - Scalar, L2-like discrete error for the flux
%
% See also convergenceUnsatBiot.

%{
Copyright 2018-2020, University of Bergen.

This file is part of the fv-unsat module.

fv-unsat is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

fv-unsat is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this file.  If not, see <http://www.gnu.org/licenses/>.
%} 


% Importing modules
mrstModule add fvbiot fv-unsat

% Creating the grid
nx = cells; % Cells in x-direction
ny = cells; % Cells in y-direction
Lx = 1; % Lenght in x-direction
Ly = 1; % Length in y-direction
G = cartGrid([nx, ny], [Lx, Ly]); % Create Cartesian grid
G = computeGeometry(G); % compute geometry

% Extracting grid information
V = G.cells.volumes; % Cell volumes
A = G.faces.areas; % Face areas
xc = G.cells.centroids(:, 1); % cell centers in x-direction
yc = G.cells.centroids(:, 2); % cell centers in y-direction
xf = G.faces.centroids(:, 1); % face centers in x-direction
yf = G.faces.centroids(:, 2); % face centers in y-direction
nu_x = G.faces.normals(:, 1); % face normals in x-direction
nu_y = G.faces.normals(:, 2); % face normals in y-direction

% Physical parameters
phys = struct(); % initializing structure to store parameters
phys.flow.mu = 1; % viscosity
phys.flow.rho = 1; % density
phys.flow.g = 1; % gravity
phys.flow.gamma = phys.flow.g * phys.flow.rho; % specific gravity
phys.flow.perm = ones(G.cells.num, 1); % permeability
phys.flow.alpha = 0.04; % vGM parameter
phys.flow.n = 2; % vGM parameter
phys.flow.m = 1 - (1 / phys.flow.n); % vGM parameter
phys.flow.theta_s = 0.4; % water content at saturation conditions
phys.flow.theta_r = 0.1; % residual water content

% Boundary and Initial Conditions

% Boundary indices
x_min = find(xf == 0); % west faces
x_max = find(xf > 0.9999*Lx & xf < 1.0001*Lx ); % east faces

y_min = find(yf == 0); % south faces 
y_max = find(yf > 0.9999*Ly & yf < 1.0001*Ly ); % north faces 

% Boundary cond. structure
boundPsi = -1;
bc = addBC([], x_min, 'pressure', boundPsi);  
bc = addBC(bc, x_max, 'pressure', boundPsi);  
bc = addBC(bc, y_min, 'pressure', boundPsi);  
bc = addBC(bc, y_max, 'pressure', boundPsi);  

% Boundary cond. values 
bcVals = zeros(G.faces.num, 1);   % initializing
bcVals(x_min) = boundPsi; % west faces
bcVals(x_max) = boundPsi; % east faces
bcVals(y_min) = boundPsi; % south faces
bcVals(y_max) = boundPsi; % north faces

% Initial condition
psi_init = boundPsi * ones(G.cells.num, 1);       

% Calling MPFA routine and creating operators
mpfa_discr = mpfa(G, phys.flow, [], 'bc', bc, 'invertBlocks', 'matlab');
     
% Time parameters
time_param = struct();  % initializing structure to store parameters
time_param.initial = 0; % initial simulation time
time_param.simTime = 1; % final simulation time  
time_param.tau = time_param.simTime/timeLevels; % constant time step
time_param.time = 0; % current time

% Retrieving analytical forms: Note that these are stored as function
% handles and retrieved from data/exactFormsRE.mat
pth = fullfile(mrstPath('fv-unsat'), 'data', 'exactFormsRE.mat');
load(pth, 'exactRE');
psi_ex = exactRE.psi;
f_ex = exactRE.source;
q_ex = exactRE.velocity;

% Calling Richards' equation model
modelEqs = modelRE(G, phys, mpfa_discr, bc, bcVals, 'arithmetic', 'off');
                           
% Time loop
solver_param = struct(); % Initializing structure to store parameters
solver_param.tol = 1E-8; % tolerance
solver_param.maxIter = 10; % maximum number of iterations

psi = psi_init; % current pressure head
time_param.time = time_param.tau; % current time

while time_param.time < time_param.simTime

    psi_n = psi; % current time level (n-index)    
    source = f_ex(time_param.time, xc, yc); % computing source term

    % Calling newton solver
    [psi, psi_m, ~] = solverRE(psi_n, modelEqs, time_param, ...
        solver_param, source);
    
    time_param.time  = time_param.time + time_param.tau; % increase time

end

% Collecting results and computing error

% True solution
psi_true = psi_ex(time_param.simTime, xc, yc); % exact pressure head
q_true = q_ex(time_param.simTime, xf, yf); % Exact velocity field
Q_true = q_true(1:G.faces.num) .* nu_x ... % Exact (normal) fluxes
         + q_true(G.faces.num+1:end) .* nu_y;

% Numerical solution
psi_num = psi; % Numerical pressure head
Q_num = modelEqs.Q(psi, psi_m); % Numerical fluxes

% Computing errors
errorPsi = sqrt(sum(V .* (psi_true - psi_num).^2)) ... 
            ./ sqrt(sum(V .* psi_true.^2));
errorFlux = sqrt(sum(A .* (Q_true - Q_num).^2)) ... 
            ./ sqrt(sum(A .* Q_true.^2));
