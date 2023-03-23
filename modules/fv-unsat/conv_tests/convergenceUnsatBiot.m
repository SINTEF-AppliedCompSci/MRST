function [e_p, e_u, e_Q, e_T] = convergenceUnsatBiot(cells, timeLevels)
% Converge of unsaturated poroelastic equations (Ex. 2 from the chapter)
%
% SYNOPSIS:
%   function [e_p, e_u, e_Q, e_T] = convergenceUnsatBiot(cells, timeLevels)
%
% PARAMETERS:
%   cells       - Scalar, number of cells = nx = ny
%   timeLevels  - Scalar, number of discrete time levels
%
% RETURNS:
%   e_p         - Scalar, L2-like discrete error for the pressure
%   e_u         - Scalar, L2-like discrete error for the displacement
%   e_Q         - Scalar, L2-like discrete error for the flux
%   e_T         - Scalar, L2-like discrete error for the traction
%
% See also convergenceRE.

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

% Creating structured triangular grid
[xx, yy] = meshgrid(0:Lx/nx:Lx, 0:Ly/ny:Ly);
tri = delaunay(xx(:), yy(:));
G = triangleGrid([xx(:) yy(:)], tri);
G = computeGeometry(G);

% Extracting grid information
V = G.cells.volumes; % Cell volumes
A = G.faces.areas; % Face areas
xc = G.cells.centroids(:, 1); % cell centers in x-direction
yc = G.cells.centroids(:, 2); % cell centers in y-direction
xf = G.faces.centroids(:, 1); % face centers in x-direction
yf = G.faces.centroids(:, 2); % face centers in y-direction
nu_x = G.faces.normals(:, 1); % face normals in x-direction
nu_y = G.faces.normals(:, 2); % face normals in y-direction
zetaf = Ly - yf; % face centers of elevation head

% Physical parameters
phys = struct(); % initialize structure to store parameters

% Flow parameters
phys.flow.rho = 1; % fluid density
phys.flow.g   = 1; % gravity acceleration
phys.flow.gamma = phys.flow.rho * phys.flow.g; % specific gravity
phys.flow.C_w   = 1; % fluid compressibility
phys.flow.mu  = 1; % fluid viscosity
phys.flow.perm = ones(G.cells.num, 1); % permeability
phys.flow.poro = 0.4; % reference porosity
phys.flow.alpha = 0.04; % vGM parameter
phys.flow.n = 2; % vGM parameter
phys.flow.m = 1 - (1 / phys.flow.n); % vGM parameter
phys.flow.theta_s = phys.flow.poro; % Water content at saturation conditions
phys.flow.theta_r = 0.1; % Residual water content
phys.flow.S_r = phys.flow.theta_r / phys.flow.poro; % Residual saturation
phys.flow.a = phys.flow.alpha / phys.flow.gamma; % vGM artificial parameter

% Mechanical parameters
phys.mech.lambda = ones(G.cells.num, 1); % First Lame paramter
phys.mech.mu = ones(G.cells.num, 1); % Second Lame parameter
phys.mech.C_s = 0.1; % solid compressibility
phys.mech.stiff = shear_normal_stress(G.cells.num, G.griddim, ...
    phys.mech.mu, phys.mech.lambda, 0 .* phys.mech.mu); % stiffness matrix

phys.flow.C_m = 1; % Porous medium compressibility
phys.flow.alpha_biot = 1 - phys.mech.C_s/phys.flow.C_m; % Biot coefficient

% Boundary and Initial Conditions

% Boundary indices
x_min = find(xf == 0); % west faces
x_max = find(xf > 0.9999*Lx & xf < 1.0001*Lx ); % east faces

y_min = find(yf == 0); % south faces 
y_max = find(yf > 0.9999*Ly & yf < 1.0001*Ly ); % north faces 

% MECHANICS BC SETUP

% Mechanics boundary cond. structure
bcMech = addBC([], x_min, 'pressure', 0);       
bcMech = addBC(bcMech, x_max, 'pressure', 0);  
bcMech = addBC(bcMech, y_min, 'pressure', 0);  
bcMech = addBC(bcMech, y_max, 'pressure', 0);      

% Mechanics boundary condition values
bcMechVals = zeros(G.griddim * G.faces.num, 1);

% FLOW BC SETUP

% Flow boundary cond. structure
boundP  = -1;
bcFlow = addBC([], x_min, 'pressure', boundP);       
bcFlow = addBC(bcFlow, x_max, 'pressure', boundP);  
bcFlow = addBC(bcFlow, y_min, 'pressure', boundP);  
bcFlow = addBC(bcFlow, y_max, 'pressure', boundP);   

% Flow boundary cond. values 
bcFlowVals = zeros(G.faces.num, 1);   
bcFlowVals(x_min) = boundP + phys.flow.gamma * zetaf(x_min);
bcFlowVals(x_max) = boundP + phys.flow.gamma * zetaf(x_max);  
bcFlowVals(y_min) = boundP + phys.flow.gamma * zetaf(y_min);  
bcFlowVals(y_max) = boundP + phys.flow.gamma * zetaf(y_max);   

% Initial conditions
p_init = -1 * ones(G.cells.num, 1);       
u_init = zeros(G.cells.num * G.griddim, 1);

% Calling MPSA/MPFA routines and creating operators

% Discretise mechanics problem
mpsa_discr = mpsa(G, phys.mech.stiff, [], 'invertBlocks', 'matlab', ...
    'bc', bcMech);

% Discretise flow problem
mpfa_discr = mpfa(G, phys.flow, [], 'invertBlocks', 'matlab', ...
    'bc', bcFlow);
     
% Time parameters
time_param = struct(); % initializing structure to store parameters
time_param.initial = 0; % initial simulation time
time_param.simTime = 1; % final simulation time  
time_param.tau = time_param.simTime/timeLevels; % constant time step
time_param.time = 0; % current time

% Retrieving analytical forms: Note that these are stored as function
% handles and retrieved from data/exactFormsREBiot.mat
pth = fullfile(mrstPath('fv-unsat'), 'data', 'exactFormsUnsatBiot.mat');
load(pth, 'exactUnsatBiot');
p_ex     = exactUnsatBiot.pressure;
fflow_ex = exactUnsatBiot.sourceFlow;
q_ex     = exactUnsatBiot.velocity;
u_ex     = exactUnsatBiot.displacement;
fmech_ex = exactUnsatBiot.sourceMech;
sxx_ex   = exactUnsatBiot.stress_xx;
syy_ex   = exactUnsatBiot.stress_yy;
sxy_ex   = exactUnsatBiot.stress_xy;

% Calling unsatBiot model
modelEqs = modelUnsatBiot(G, phys, mpfa_discr, mpsa_discr,...
    bcFlow, bcFlowVals, bcMech, bcMechVals, 'arithmetic', 'on');

% Time loop
solver_param = struct(); % Initializing structure to store parameters
solver_param.tol = 1E-8; % tolerance
solver_param.maxIter = 10; % maximum number of iterations

p = p_init; % current pressure
u = u_init; % current displacement
time_param.time = time_param.tau; % current time

while time_param.time < time_param.simTime
    
    p_n = p; % current time level (n-index) 
    u_n = u; % current time level (n-index)
    
    % Source terms
    sourceMech = zeros(G.cells.num * G.griddim, 1);
    f_mech = fmech_ex(time_param.time, xc, yc);
    sourceMech(1:G.griddim:end) = f_mech(1:G.cells.num);
    sourceMech(2:G.griddim:end) = f_mech(G.cells.num+1:end);
    sourceFlow = fflow_ex(time_param.time, xc, yc);

    % Calling Newton solver
    [p, p_m, u, ~] = solverUnsatBiot(G, p_n, u_n, modelEqs, time_param, ...
        solver_param, sourceFlow, sourceMech);
    
    time_param.time  = time_param.time + time_param.tau; % increase time

end

% Collecting results and computing errors

% Exact pressure
p_true = p_ex(time_param.simTime, xc, yc);              

% Exact displacement field
u_true = zeros(G.cells.num * G.griddim, 1);   
displacement = u_ex(time_param.simTime, xc, yc);               
u_true(1:G.griddim:end) = displacement(1:G.cells.num);
u_true(2:G.griddim:end) = displacement(G.cells.num+1:end);

% Exact normal fluxes
q_true = q_ex(time_param.simTime, xf, yf);               
Q_true = q_true(1:G.faces.num).*nu_x + q_true(G.faces.num+1:end).*nu_y; 

% Exact traction
sxx = sxx_ex(time_param.simTime, xf, yf);
syy = syy_ex(time_param.simTime, xf, yf);
sxy = sxy_ex(time_param.simTime, xf, yf);
T_true = zeros(G.faces.num * G.griddim, 1);
T_true(1:G.griddim:end) = sxx .* nu_x + sxy .* nu_y;
T_true(2:G.griddim:end) = sxy .* nu_x + syy .* nu_y;

% Numerical solutions
p_num = p; % Numerical pressure 
u_num = u; % Numerical displacement
Q_num = modelEqs.Q(p, p_m); % Numerical fluxes
T_num = modelEqs.T(u); % Numerical traction    

% Computing errors
e_p = sqrt(sum(V .* (p_true - p_num).^2)) ...
    ./ sqrt(sum(V .* p_true.^2));
e_u = sqrt(sum(modelEqs.sca2vec(V) .* (u_true - u_num).^2)) ...
    ./ sqrt(sum(modelEqs.sca2vec(V) .* u_true.^2));             
e_Q = sqrt(sum(A .* (Q_true - Q_num).^2)) ...
    ./ sqrt(sum(A .* Q_true.^2));
e_T = sqrt(sum(modelEqs.sca2vec(A) .* (T_true - T_num).^2)) ...
    ./ sqrt(sum(modelEqs.sca2vec(A) .* T_true.^2)); 