%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               fv-unsat                                                
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Example 3: Water inflitration in an initially dry soil
% Author: Jhabriel Varela. E-mail: Jhabriel.Varela@uib.no. 

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


%% Importing required modules
clear; clc(); mrstModule add fvbiot

%% Setting up the grid
nx = 5;     ny = 5;     nz = 30;     % cells        
Lx = 1;     Ly = 1;     Lz = 1;      % domain lenght [m]       
G = cartGrid([nx,ny,nz],[Lx,Ly,Lz]); % create Cartesian Grid
G = computeGeometry(G);              % compute geometry

%Plotting grid
newplot; plotGrid(G); axis off; 
pbaspect([1,1,5]); view([-51,26]);

%% Physical properties 
soil = getHydraulicProperties('newMexSample');
rho = 1 * gram / (centi * meter)^3;          % [kg/m^3] water density
mu = 0.01 * gram / (centi * meter * second); % [kg/(m*s)] water viscosity
g = 980.66 * centi * meter / (second^2);     % [kg/(m*s^2)] gravity acceleration
gamma = rho*g;                               % [kg/(m^2*s^2)] specific gravity
K_sat = soil.K_s;                            % [m/s] saturated hydraulic conductivity
k = (K_sat*mu)/(rho*g);                      % [m^2] intrinsic permeability
rock.perm = k * ones(G.cells.num, 1);        % creating perm structure
alpha = soil.alpha;                          % [1/m] Equation parameter
nVan = soil.n;                               % [-] Equation parameter
mVan = 1-(1/nVan);                           % [-] Equation parameter
theta_s = soil.theta_s;                      % [-] Saturation soil moisture
theta_r = soil.theta_r;                      % [-] Residual soil moisture

%% Water retention curves
[theta, krw, C_theta] = vanGenuchtenMualemTheta(alpha, theta_s, ...
                                                theta_r, nVan, mVan);

%% Boundary and Initial Conditions

% Extracting grid information
zc = G.cells.centroids(:, 3);  % cell centers in z-direction
zf = G.faces.centroids(:, 3);  % face centers in z-direction
zetac = Lz - zc;               % centroids of cells of elev. head
zetaf = Lz - zf;               % centroids of faces of elev. head
V = G.cells.volumes;           % Cell volumes
z_min = find(zf == 0);         % top faces idx
z_max = find(zf > 0.9999*Lz & zf < 1.0001*Lz );  % bottom faces idx

% Creating the boundary structure
psiT = -75 * centi * meter;    % [m] Top boundary pressure head
psiB = -1000 * centi * meter;  % [m] Bottom boundary pressure head
bc = addBC([], z_min, 'pressure', psiT);       
bc = addBC(bc, z_max, 'pressure', psiB);      
bcVal = zeros(G.faces.num, 1);
bcVal(z_min) = psiT + zetaf(z_min); % assigning Top boundary
bcVal(z_max) = psiB + zetaf(z_max); % assigning Bottom boundary
% Note that we are adding the gravity contributions to bcVal. This must
% be done for every Dirichlet face only to bcVal, not to bc.

% Initial Condition
psi_init = psiB .* ones(G.cells.num, 1);  % [m] Initial condition

%% Calling MPFA routine and creating discrete operators
mpfa_discr = mpfa(G, rock, [], 'bc', bc, 'invertBlocks', 'matlab');

% Flux discretization
F       = @(x) mpfa_discr.F * x;         
% Boundary discretization
boundF  = @(x) mpfa_discr.boundFlux * x; 
% Discrete divergence
divF    = @(x) mpfa_discr.div * x;       

%% Creating AD variable
psi_ad = initVariablesADI(psi_init);  % Initial value set to psi_init

%% Time parameters
simTime = 72 * hour; % [s] final simulation time
tau_init = 10;       % [s] initial time step
tau_min = 10;        % [s] minimum time step 
tau_max = 10000;     % [s] maximum time step
tau = tau_init;      % [s] initializing time step
timeCum = 0;         % [s] initializing cumulative time

%% Printing parameters
printLevels = 20;       % number of printed levels
printTimes = ((simTime/printLevels):(simTime/printLevels):simTime)';
pp = 1;                 % initializing print counter
ee = 1;                 % initializing export counter

%% Discrete equations

% Arithmetic average of krw
krwAr = @(psi_m) arithmeticAverageMPFA(G, bc, krw, psi_m);

% Darcy Flux
Q = @(psi, psi_m) (rho.*g./mu) .* krwAr(psi_m) .* ...
    (F(psi + zetac) + boundF(bcVal));

% Mass Conservation                         
psiEq = @(psi, psi_n, psi_m, tau, source)   (V./tau) .* (theta(psi_m) ...
    + C_theta(psi_m) .* (psi - psi_m) - theta(psi_n)) ...
    + divF(Q(psi, psi_m)) - V .* source;
                            
%% Creating solution structure
sol.time  = zeros(printLevels,1); % time
sol.psi   = cell(printLevels,1);  % pressure head
sol.theta = cell(printLevels,1);  % water content
sol.flux  = cell(printLevels,1);  % flux
                               
%% Time loop
tol = 1E-6; % tolerance
maxIter = 10; % maximum number of iterations
while timeCum < simTime
    
    psi_n = psi_ad.val;             % [m] current time step h (n-index)
    timeCum = timeCum + tau;        % [s] cumulative time
    source = zeros(G.cells.num,1);  % source term equal to zero
            
    % Newton loop
    [psi_ad, psi_m, iter] = solverRE(psi_ad, psi_n, psiEq, tau, ...
        source, timeCum, tol, maxIter);
                    
    % Time stepping routine
    [tau, pp] = timeStepping(tau, tau_min, tau_max, simTime, timeCum, ...
        iter, printTimes, pp);
    
    % Storing solutions at each printing time
    if timeCum == printTimes(ee)
        sol.time(ee,1) = timeCum;
        sol.psi{ee,1} = psi_ad.val;
        sol.theta{ee,1} = theta(psi_ad.val);
        sol.flux{ee,1} = Q(psi_ad.val, psi_m);
        ee = ee + 1;   
    end
    
end

fprintf('\n sol: \n'); disp(sol); % printing sol structure in console

%% Plotting pressure head and water content profiles 
for ii=1:printLevels
    newplot;
    subplot(1,2,1); % pressure head plots
    plotCellData(G, sol.psi{ii,1} / centi);
    title('\psi_w [cm]', 'FontSize', 13)
    axis off; pbaspect([1,1,5]); view([-51,26]); 
    cb = colorbar; cb.FontSize=11;
    t_h1=handle(text('Units','normalized',...     % Use [%] of the axis length
                     'Position',[0.5,-0.05,1],... % Position of text   
                     'EdgeColor','w'));           % Textbox 
    t_h1.String=sprintf('Time: %2.1f [hours]',sol.time(ii)/hour); 
    set(gca,'Ydir','reverse'); axis tight; box on; grid on; 
    
    subplot(1,2,2); % water content plots
    plotCellData(G, sol.theta{ii,1})
    title('\theta_w [ - ]', 'FontSize', 13)
    axis off; pbaspect([1,1,5]); view([-51,26]);
    cb = colorbar; cb.FontSize=11;
    t_h2=handle(text('Units','normalized',...     % Use [%] of the axis length
                     'Position',[0.5,-0.05,1],... % Position of text   
                     'EdgeColor','w'));           % Textbox 
    t_h2.String=sprintf('Time: %2.1f [hours]',sol.time(ii)/hour);
    set(gca,'Ydir','reverse'); axis tight; box on; grid on; 
    
    pause(0.01); t_h1.String =[]; t_h2.String =[];
end