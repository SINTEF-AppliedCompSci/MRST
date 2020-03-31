%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               fv-unsat                                  %              
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
clear; clc(); mrstModule add fvbiot fv-unsat

%% Setting up the grid
nx = 5;     ny = 5;     nz = 30;          % cells        
Lx = 1;     Ly = 1;     Lz = 1;           % domain lenght [m]       
G = cartGrid([nx, ny, nz], [Lx, Ly, Lz]); % create Cartesian Grid
G = computeGeometry(G);                   % compute geometry

%Plotting grid
newplot; plotGrid(G); axis off; 
pbaspect([1, 1, 5]); view([-51, 26]);

%% Physical properties 
soil = getHydraulicProperties('newMexSample');  % get soil properties
phys = struct(); % create structure to store physical properties

% Flow parameters
phys.flow.rho = 1 * gram / (centi * meter)^3; % density
phys.flow.mu = 0.01 * gram / (centi * meter * second); % viscosity
phys.flow.g = 980.66 * centi * meter / (second^2); % gravity 
phys.flow.gamma = phys.flow.rho * phys.flow.g; % specific gravity
phys.flow.K = soil.K_s; % saturated hydraulic conductivity
phys.flow.perm = (phys.flow.K * phys.flow.mu / phys.flow.gamma) .* ...
    ones(G.cells.num, 1); % intrinsic permeability
phys.flow.alpha = soil.alpha / meter; % vGM parameter
phys.flow.n = soil.n; % vGM parameter
phys.flow.m = 1-(1/phys.flow.n); % vGM parameter
phys.flow.theta_s = soil.theta_s; % Water content at saturation conditions
phys.flow.theta_r = soil.theta_r; % Residual water content

%% Boundary and Initial Conditions

% Extracting grid information
zf = G.faces.centroids(:, end); % face centers in z-direction
zetaf = Lz - zf; % centroids of faces of elev. head
z_min = find(zf == 0); % top faces idx
z_max = find(zf > 0.9999*Lz & zf < 1.0001*Lz ); % bottom faces idx

% Creating the boundary structure
psiT = -75 * centi * meter; % Top boundary pressure head
psiB = -1000 * centi * meter; % Bottom boundary pressure head
bc = addBC([], z_min, 'pressure', psiT);       
bc = addBC(bc, z_max, 'pressure', psiB);      
bcVal = zeros(G.faces.num, 1);
bcVal(z_min) = psiT + zetaf(z_min); % assigning Top boundary
bcVal(z_max) = psiB + zetaf(z_max); % assigning Bottom boundary
% Note that we are adding the gravity contributions to bcVal. This must
% be done for every Dirichlet face only to bcVal, not to bc.

% Initial Condition
psi = psiB .* ones(G.cells.num, 1);

%% Discretize the flow problem using MPFA
mpfa_discr = mpfa(G, phys.flow, [], 'bc', bc, 'invertBlocks', 'matlab');

%% Time parameters
time_param = struct(); % Initialize structure to store time params
time_param.simTime = 72 * hour; % final simulation time
time_param.tau_init = 10 * second; % initial time step
time_param.tau_min = 10 * second; % minimum time step 
time_param.tau_max = 10000 * second; % maximum time step
time_param.tau = time_param.tau_init; % initializing time step
time_param.time = 0; % initializing current time

%% Printing parameters
print_param = struct(); % Initialize structure to store time params
print_param.levels = 20; % number of printed levels
print_param.times  = ((time_param.simTime/print_param.levels) : ...
    (time_param.simTime/print_param.levels):time_param.simTime)';
print_param.print = 1; % initializing print counter
print_param.export = 1; % intializing export counter

%% Call Richards' equation model 
modelEqs = modelRE(G, phys, mpfa_discr, bc, bcVal, 'arithmetic', 'on'); 
                         
%% Creating solution structure
sol.time = zeros(print_param.levels, 1);  
sol.psi = cell(print_param.levels, 1);   
sol.theta = cell(print_param.levels, 1);   
sol.flux = cell(print_param.levels, 1);  
                               
%% Time loop
solver_param = struct(); % initialize structure to store solver parameters
solver_param.tol = 1E-6; % tolerance
solver_param.maxIter = 10; % maximum number of iterations

while time_param.time < time_param.simTime
    
    psi_n = psi; % current time step (n-index)
    time_param.time = time_param.time + time_param.tau; % current time
    source = zeros(G.cells.num,1); % source term equal to zero
            
    % Newton loop
    [psi, psi_m, iter] = solverRE(psi_n, modelEqs, time_param, ...
        solver_param, source);
                    
    % Determine next time step
    [time_param.tau, print_param.print] = timeStepping(time_param, ...
        print_param, iter);
    
    % Storing solutions at each printing time
    if time_param.time == print_param.times(print_param.export)
        sol.time(print_param.export,1) = time_param.time;
        sol.psi{print_param.export,1} = psi;
        sol.theta{print_param.export,1} = modelEqs.theta(psi);
        sol.flux{print_param.export,1} = modelEqs.Q(psi, psi_m);
        print_param.export = print_param.export + 1;   
    end
    
end

% Displaying "sol" structure in console
fprintf('\n sol: \n'); disp(sol); 

%% Plotting pressure head and water content profiles 
for ii=1:print_param.levels
    newplot;
    subplot(1,2,1); % Pressure head plots
    plotCellData(G, sol.psi{ii,1} / centi);
    title('\psi_w [cm]', 'FontSize', 13)
    axis off; pbaspect([1,1,5]); view([-51,26]); 
    cb = colorbar; cb.FontSize=11;
    t_h1=handle(text('Units','normalized', ... % Use [%] of the axis length
                     'Position',[0.5,-0.05,1], ... % Position of text   
                     'EdgeColor','w')); % Textbox 
    t_h1.String=sprintf('Time: %2.1f [hours]', sol.time(ii) / hour); 
    set(gca,'Ydir','reverse'); axis tight; box on; grid on; 
    
    subplot(1,2,2); % Water content plots
    plotCellData(G, sol.theta{ii,1})
    title('\theta_w [ - ]', 'FontSize', 13)
    axis off; pbaspect([1,1,5]); view([-51,26]);
    cb = colorbar; cb.FontSize=11;
    t_h2=handle(text('Units','normalized', ... % Use [%] of the axis length
                     'Position',[0.5,-0.05,1], ... % Position of text   
                     'EdgeColor','w')); % Textbox 
    t_h2.String=sprintf('Time: %2.1f [hours]', sol.time(ii) / hour);
    set(gca,'Ydir','reverse'); axis tight; box on; grid on; 
    
    pause(0.01); t_h1.String =[]; t_h2.String =[];
end