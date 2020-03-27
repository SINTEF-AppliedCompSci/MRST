function [e_p, e_u, e_Q, e_T] = convergence_re_biot(cells, timeLevels)
% Converge of unsaturated poroelastic equations (Ex. 2 from the chapter)
%
% SYNOPSIS:
%function [error_p, error_u, error_flux, error_traction] ...
%    = convergence_re_biot(cells, timeLevels)
%
% PARAMETERS:
%   cells       - Scalar, number of cells = nx = ny.
%   timeLevels  - Scalar, number of discrete time levels.
%
%  RETURNS:
%   e_p         - Scalar, L2-like discrete error for the pressure.
%   e_u         - Scalar, L2-like discrete error for the displacement
%   e_Q         - Scalar, L2-like discrete error for the flux.
%   e_T         - Scalar, L2-like discrete error for the traction.
%

%{
Copyright 2018-2019, University of Bergen.

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
mrstModule add fvbiot

% Creating the grid
nx = cells;  % Cells in x-direction
ny = cells;  % Cells in y-direction
Lx = 1;      % Lenght in x-direction
Ly = 1;      % Length in y-direction

% Creating structured triangular grid
[xx, yy] = meshgrid(0:Lx/nx:Lx,0:Ly/ny:Ly);
tri      = delaunay(xx(:), yy(:));
G        = triangleGrid([xx(:) yy(:)], tri);
G        = computeGeometry(G);

% Extracting grid information
V = G.cells.volumes;           % Cell volumes
A = G.faces.areas;             % Face areas
xc = G.cells.centroids(:, 1);  % cell centers in x-direction
yc = G.cells.centroids(:, 2);  % cell centers in y-direction
xf = G.faces.centroids(:, 1);  % face centers in x-direction
yf = G.faces.centroids(:, 2);  % face centers in y-direction
nu_x = G.faces.normals(:, 1);  % face normals in x-direction
nu_y = G.faces.normals(:, 2);  % face normals in y-direction
zetac = Ly - yc;               % cell centers of elevation head
zetaf = Ly - yf;               % face centers of elevation head

% Physical parameters

% Flow parameters
rho_w = 1;              % fluid density
g = 1;                  % gravity
gamma = rho_w * g;      % specific gravity
C_w = 1;                % fluid compressibility
mu_w = 1;               % fluid dynamics viscosity
k = 1;                  % intrinsic permeability
rock.perm = k * ones(G.cells.num, 1); % creating perm structure
n = 0.4;                % reference porosity

% Water retention curves (van Genuchten - Mualem)
alpha = 0.04;           % Equation parameter
nVan = 2;               % Equation parameter
mVan = 1-(1/nVan);      % Equation parameter
theta_r = 0.1;          % Residual water content 
S_r = theta_r/n;        % Residual water saturation
a_v = alpha / gamma;    % Artifical parameter 

% Mechanical parameters
lambda_s = 1;           % first Lame parameter
mu_s = 1;               % second Lame parameter
C_s = 0.1;              % solid compressibility
constit = shear_normal_stress(G.cells.num, G.griddim, ... %
    ones(G.cells.num, 1) .* mu_s, ...                     % creating stiff
    ones(G.cells.num, 1) .* lambda_s, ...                 % matrix
    zeros(G.cells.num, 1));                               %

% Coupling parameters
C_m = 1;                  % porous medium compressibility
alpha_biot = 1 - C_s/C_m; % Biot's coefficient

% Water retention curves

[S_w, krw, C_S] = vanGenuchtenMualemSw(a_v, S_r, nVan, mVan);

% Boundary and Initial Conditions

% Boundary indices
x_min = find(xf == 0);                           % west faces
x_max = find(xf > 0.9999*Lx & xf < 1.0001*Lx );  % east faces

y_min = find(yf == 0);                           % south faces 
y_max = find(yf > 0.9999*Ly & yf < 1.0001*Ly );  % north faces 

% MECHANICS BC SETUP

% Mechanics boundary cond. structure
bcMech = addBC([], x_min, 'pressure', 0);       
bcMech = addBC(bcMech, x_max, 'pressure', 0);  
bcMech = addBC(bcMech, y_min, 'pressure', 0);  
bcMech = addBC(bcMech, y_max, 'pressure', 0);      

% Mechanics boundary condition values
bcMechVals = zeros(G.griddim * G.faces.num, 1);
bcMechVals(2*x_min-1) = 0; % ux = 0 at west bound
bcMechVals(2*x_min)   = 0; % uy = 0 at west bound
bcMechVals(2*x_max-1) = 0; % ux = 0 at east bound
bcMechVals(2*x_max)   = 0; % uy = 0 at east bound
bcMechVals(2*y_min-1) = 0; % ux = 0 at south bound
bcMechVals(2*y_min)   = 0; % uy = 0 at south bound
bcMechVals(2*y_max-1) = 0; % ux = 0 at north bound
bcMechVals(2*y_max)   = 0; % uy = 0 at north bound

% FLOW BC SETUP

% Flow boundary cond. structure
boundP  = -1;
bcFlow = addBC([], x_min, 'pressure', boundP);       
bcFlow = addBC(bcFlow, x_max, 'pressure', boundP);  
bcFlow = addBC(bcFlow, y_min, 'pressure', boundP);  
bcFlow = addBC(bcFlow, y_max, 'pressure', boundP);   

% Flow boundary cond. values 
bcFlowVals = zeros(G.faces.num, 1);   
bcFlowVals(x_min) = boundP + gamma * zetaf(x_min);
bcFlowVals(x_max) = boundP + gamma * zetaf(x_max);  
bcFlowVals(y_min) = boundP + gamma * zetaf(y_min);  
bcFlowVals(y_max) = boundP + gamma * zetaf(y_max);   

% Initial conditions
p_init = -1 * ones(G.cells.num, 1);       
u_init = zeros(G.cells.num * G.griddim, 1);

% Calling MPSA/MPFA routines and creating operators

% Discretise mechanics problem
mpsa_discr = mpsa(G, constit, [], 'invertBlocks', 'matlab', 'bc', bcMech);

% Discretise flow problem
mpfa_discr = mpfa(G, rock, [], 'invertBlocks', 'matlab', 'bc', bcFlow);
     
% MPSA discrete operators
S = @(x) mpsa_discr.stress * x;             % stress discretization 
boundS = @(x) mpsa_discr.boundStress * x;   % boundary mech. discretization
gradP = @(x) mpsa_discr.gradP * x;          % gradient of pressure
divU = @(x) mpsa_discr.divD * x;            % divergence of displacement
divS = @(x) mpsa_discr.div * x;             % divergence of stress
compat = @(x) mpsa_discr.stabDelta * x;  % stability parameter

% MPFA discrete operators
F = @(x) mpfa_discr.F * x;                  % flux discretization
boundF = @(x) mpfa_discr.boundFlux * x;     % boundary flow discretization
divF = @(x) mpfa_discr.div * x;             % divergence of the flux

% Time parameters
iniTime = 0;  % intial simulation time
simTime = 1;  % final simulation time 
times = linspace(iniTime, simTime, timeLevels + 1); % evaluation times
tau = diff(times); % time steps

% Retrieving analytical forms: Note that these are stored as function
% handles and retrieved from data/exactFormsREBiot.mat

load('exactFormsREBiot.mat','exactREBiot')
p_ex     = exactREBiot.pressure;
fflow_ex = exactREBiot.sourceFlow;
q_ex     = exactREBiot.velocity;
u_ex     = exactREBiot.displacement;
fmech_ex = exactREBiot.sourceMech;
sxx_ex   = exactREBiot.stress_xx;
syy_ex   = exactREBiot.stress_yy;
sxy_ex   = exactREBiot.stress_xy;

% Discrete equations

% Function that maps scalar cell-center variables to vector cell-centers
sca2vec = @(scalar) repmat(scalar, [G.griddim, 1]); 

% Mechanics equations
T    = @(u)               S(u) + boundS(bcMechVals);
uEq1 = @(u)               divS(T(u));
uEq2 = @(p, p_n, f_mech)  -alpha_biot .* gradP(p.*S_w(p_n)) ...
                          +sca2vec(V) .* f_mech;
% Flow equations
krwAr = @(p_m) arithmeticAverageMPFA(G, bcFlow, krw, p_m);
Q     = @(p,p_m) (krwAr(p_m) ./ mu_w) .* (F(p + gamma .* zetac) + ...
        boundF(bcFlowVals));
xi    = @(p_n) (alpha_biot-n).*C_s.*S_w(p_n).^2 + n.*C_w.*S_w(p_n);
chi   = @(p_n) (alpha_biot-n).*C_s.*S_w(p_n).*p_n + n;  
pEq1  = @(p_n,u,u_n) alpha_biot .* S_w(p_n) .* divU(u-u_n);
pEq2  = @(p, p_n, p_m, tau, f_flow)  ...
        alpha_biot.^2 .* S_w(p_n) .* compat(S_w(p_n) .* p_n) ...
        + V .* xi(p_n) .* (p - p_n)  ...
        + V .* chi(p_n) .* (S_w(p_m) + C_S(p_m).*(p - p_m) - S_w(p_n)) ...
        + tau .* divF(Q(p, p_m)) - V .* tau .* f_flow;

% Time loop
tt = 1; % time counter
tol = 1E-8; % tolerance
maxIter = 10; % maximum number of iterations
p = p_init; % current pressure
u = u_init; % current displacement
while times(tt) < simTime
    
    p_n = p;      % current time level (n-index) 
    u_n = u;      % current time level (n-index)
    tt = tt + 1;  % increasing time counter
    
    % Source terms
    sourceMech = zeros(G.cells.num * G.griddim, 1);
    f_mech = fmech_ex(times(tt), xc, yc);
    sourceMech(1:G.griddim:end) = f_mech(1:G.cells.num);
    sourceMech(2:G.griddim:end) = f_mech(G.cells.num+1:end);
    sourceFlow = fflow_ex(times(tt), xc, yc);

    % Calling Newton solver
    [p, p_m, u, ~] = solverUnsatBiot(G, p_n, u_n, ...
        pEq1, pEq2, uEq1, uEq2, tau(tt-1), sourceFlow, sourceMech, ...
        times(tt), tol, maxIter);
    
end

% Collecting results and computing errors

% Exact pressure
p_true = p_ex(simTime, xc, yc);              

% Exact displacement field
u_true = zeros(G.cells.num * G.griddim, 1);   
displacement = u_ex(simTime, xc, yc);               
u_true(1:G.griddim:end) = displacement(1:G.cells.num);
u_true(2:G.griddim:end) = displacement(G.cells.num+1:end);

% Exact normal fluxes
q_true = q_ex(simTime, xf, yf);               
Q_true = q_true(1:G.faces.num).*nu_x + q_true(G.faces.num+1:end).*nu_y; 

% Exact traction
sxx = sxx_ex(simTime, xf, yf);
syy = syy_ex(simTime, xf, yf);
sxy = sxy_ex(simTime, xf, yf);
T_true = zeros(G.faces.num * G.griddim, 1);
T_true(1:G.griddim:end) = sxx .* nu_x + sxy .* nu_y;
T_true(2:G.griddim:end) = sxy .* nu_x + syy .* nu_y;

% Numerical solutions
p_num = p;          % Numerical pressure 
u_num = u;          % Numerical displacement
Q_num = Q(p, p_m);  % Numerical fluxes
T_num = T(u);       % Numerical traction    

% Computing errors
e_p = sqrt(sum(V .* (p_true - p_num).^2)) ./ sqrt(sum(V .* p_true.^2));
e_u = sqrt(sum(sca2vec(V) .* (u_true - u_num).^2)) ./ sqrt(sum(sca2vec(V) .* u_true.^2));             
e_Q = sqrt(sum(A .* (Q_true - Q_num).^2)) ./ sqrt(sum(A .* Q_true.^2));
e_T = sqrt(sum(sca2vec(A) .* (T_true - T_num).^2)) ./ sqrt(sum(sca2vec(A) .* T_true.^2)); 