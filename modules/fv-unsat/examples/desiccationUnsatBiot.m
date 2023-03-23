%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               fv-unsat                                  %              
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Example 4: Desiccation of a clay sample in a Petri-dish
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
clear; clc(); mrstModule add fv-unsat fvbiot distmesh

%% Creating grid

% Two-dimensional grid
r = 50 * milli * meter; % radii of the Petri-dish
fd =@(p) sqrt(sum(p.^2, 2)) - r; % circular domain function
min_x = -r;  max_x = r; % min and max values in x-axis
min_y = -r;  max_y = r; % min and max values in y-axis                  
h = (2*r)/25; % step size
[p, t] = distmesh2d(fd, @huniform, h, [min_x, min_y; max_x, max_y], []);
p = p + r; % shifting triangulation points
G = triangleGrid(p, t); % creating triangular grid
G = pebi(G); % creaing Voronoi diagram

% Extrude in the z-direction
Lz = 15 * milli * meter; % thickness of the Petri-dish
nz = 5; % number of layers in z-axis
dz = Lz/nz; % thickness of each layer
thick = dz .* ones(nz, 1); % thickness vector
G = makeLayeredGrid(G, thick); % extrude grid
G = computeGeometry(G); % compute geometry

% Plotting grid
newplot; plotGrid(G, 'FaceColor', [.8, .8, .8]); 
axis off equal; view([-9, 33]); 

%% Extracting grid information
Nc = G.cells.num; % total number of cells
Nf = G.faces.num; % total number of faces
Nd = G.griddim; % grid dimension
A = G.faces.areas; % face areas
xc = G.cells.centroids(:, 1); % cell centers in x-direction
yc = G.cells.centroids(:, 2); % cell centers in y-direction 
zf = G.faces.centroids(:, 3); % face centers in z-direction
zetaf = Lz - zf; % face centers of elev. head

fNei = G.faces.neighbors; % extracting faces neighbors
ext_f = find(((fNei(:,1) == 0) + (fNei(:,2) == 0)) == 1); % external faces
z_min = find(zf == 0); % idx of top faces
z_max = find(zf > 0.9999*Lz & zf < 1.0001*Lz );  % idx of bottom faces
sides = ext_f(~ismember(ext_f, [z_min; z_max])); % idx of sides faces

%% Physical parameters
phys = struct(); % create structure to store physical properties

% Due to the low permeability of the clay, we need to use a simple scaling
% technique to avoid singular matrices while inverting local mpfa systems
scale = 1E10;

% Mechanics parameters [Kaolinite]
phys.mech.lambda = 1.229E11 .* ones(Nc, 1) * Pascal; % first Lame parameter
phys.mech.mu = 4.7794E10 .* ones(Nc, 1) * Pascal; % second Lame parameter
phys.mech.C_s = 5.618E-11 / Pascal; % solid compressibility
phys.mech.rho = 1769 * kilo * gram / meter^3; % solid density
phys.mech.stiff = shear_normal_stress(Nc, Nd, ... % stiffnes matrix
    phys.mech.mu, phys.mech.lambda, 0 .* phys.mech.mu);    

% Flow parameters [Water]
soil = getHydraulicProperties('clay'); % getting clay hydraulic properties
phys.flow.K = soil.K_s * meter / second; % saturated hydraulic conduct.
phys.flow.rho = 1014 * kilo * gram / meter^3; % fluid density
phys.flow.g = 9.8006 * meter / second^2; % gravity acceleration
phys.flow.gamma = phys.flow.rho * phys.flow.g;  % specific gravity
phys.flow.mu = 1 * centi * poise * scale; % fluid dynamic viscosity
phys.flow.perm = (phys.flow.K * phys.flow.mu / phys.flow.gamma) .* ...
    ones(G.cells.num, 1); % intrinsic permeability
phys.flow.poro = soil.theta_s; % porosity
phys.flow.C_w = 4.5455E-10 / Pascal; % fluid compressibility
phys.flow.alpha = soil.alpha / meter; % vGM parameter
phys.flow.a = phys.flow.alpha / phys.flow.gamma; % Artificial parameter
phys.flow.n = soil.n; % vGM parameter
phys.flow.m = 1 - ( 1 / phys.flow.n); % vGM parameter
phys.flow.theta_r = soil.theta_r; % Residual water content
phys.flow.S_r = phys.flow.theta_r / phys.flow.poro; % Residual saturation
phys.flow.C_m = 2.17E-10 / Pascal; % Porous medium compressibility
phys.flow.alpha_biot = 1 - phys.mech.C_s/phys.flow.C_m; % Biot coefficient

% Obtaining critical pressure
phys.flow.temperature = 298.15 * Kelvin; % Ambient temperature
phys.flow.relativeHumidity = 0.5; % Ambient relative humidity
p_crit = computeCriticalPressure(phys); 

%% Boundary conditions

% MECHANICS BOUNDARY CONDITIONS

% Creating the boundary structure for the mechanics problem
bcMech = addBC([], sides, 'pressure', 0); % u=0 at the sides
bcMech = addBC(bcMech, z_max, 'pressure', 0); % u=0 at the bottom

% Boundary values are assigned to a vector of size Nd * Nf. Face
% number one will have its x-condition assigned to bcVals(1), y-condition
% at bcVals(2), and z-condition at bcVals(3). Next face number 2 etc.
% Face i has traction (or displacement) index 3*i-2, 3*i-1, 3*i (3D).
% Mixed boundary conditions, i.e., roller boundary conditions are not
% implemented
bcMechVals = zeros(Nd * Nf, 1);

% Setting lateral faces as zero displacement 
bcMechVals(Nd * sides-2) = 0; % ux = 0 at the sides
bcMechVals(Nd * sides-1) = 0; % uy = 0 at the sides
bcMechVals(Nd * sides) = 0;   % uz = 0 at the sides

% Setting bottom layer as zero displacement
bcMechVals(Nd * z_max-2) = 0; % ux = 0 at the bottom
bcMechVals(Nd * z_max-1) = 0; % uy = 0 at the bottom
bcMechVals(Nd * z_max) = 0;   % uz = 0 at the bottom

% FLOW BOUNDARY CONDITONS

% FLUX CONTROLLED BOUNDARY CONDITIONS

% We create a vector containing the indidual fluxes for each
% top face. Note that we mulitply by the viscosity since the 
% treatment of the boundary conditions for the mpfa routine assumes
% unit viscosity.
vMax = 6E-07 * meter / second; % maximum evaporation velocity
Qtop_f = vMax .* phys.flow.mu .* A(z_min);  % Top boundary flux

% Creating the boundary structure for flux-controlled BC
bcFlow_f = addBC([], z_min, 'flux', Qtop_f);     
bcFlowVals_f = zeros(Nf, 1); 
bcFlowVals_f(z_min) = Qtop_f;                        

% PRESSURE CONTROLLED BOUNDARY CONDITIONS

% Creating the boundary structure for pressure-controlled BC
bcFlow_p = addBC([], z_min, 'pressure', p_crit);  
bcFlowVals_p = zeros(Nf, 1);                          
bcFlowVals_p(z_min) = p_crit + phys.flow.gamma .* zetaf(z_min);                             
% The second term of the above equation represents the gravity 
% contribution, and must be included for Dirichlet boundary conditions

% INITIAL CONDITIONS
u_init = zeros(Nd * Nc, 1) * meter; 
p_init = -0.1 * kilo * Pascal * ones(Nc, 1);

%% Calling MPSA/MPFA routines and creating discrete operators

% Discretize mechanics problem
mpsa_discr = mpsa(G, phys.mech.stiff, [], 'invertBlocks', 'matlab', ...
    'bc', bcMech);

% Discretize flow problem for flux-controlled boundary conditions
mpfa_discr_flux = mpfa(G, phys.flow, [], 'invertBlocks', 'matlab', ...
    'bc', bcFlow_f);

% Discretize flow problem for pressure-controlled boundary conditions
mpfa_discr_pres = mpfa(G, phys.flow, [], 'invertBlocks', 'matlab', ...
    'bc', bcFlow_p);

%% Time parameters
time_param = struct(); % initialize structure to store parameters
time_param.simTime = 2 * hour; % final simulation time
time_param.tau_init = 25 * second; % initial time step
time_param.tau_min = 25 * second; % minimum time step
time_param.tau_max = 1000 * second; % maximum time step
time_param.tau = time_param.tau_init; % initializing time step
time_param.time = 0; % initializing cumulative time

%% Printing parameters
print_param = struct(); % initialize structure to store parameters
print_param.levels = 50; % number of printed levels
print_param.times  = ((time_param.simTime/print_param.levels) : ...
    (time_param.simTime/print_param.levels):time_param.simTime)';
print_param.print = 1; % initializing print counter
print_param.export = 1; % intializing export counter

%% Calling the model for the unsaturated poroelastic equations

% Setting up model for flux-controlled problem
modelEqsFlux = modelUnsatBiot(G, phys, mpfa_discr_flux, mpsa_discr, ...
    bcFlow_f, bcFlowVals_f, bcMech, bcMechVals, 'upstream', 'on');

% Setting up model for pressure-controlled problem
modelEqsPres = modelUnsatBiot(G, phys, mpfa_discr_pres, mpsa_discr, ...
    bcFlow_p, bcFlowVals_p, bcMech, bcMechVals, 'upstream', 'on');

%% Solution structure
sol.time = zeros(print_param.levels, 1); 
sol.u = cell(print_param.levels, 1); 
sol.p = cell(print_param.levels, 1); 
sol.pTop = zeros(print_param.levels, 1); 
sol.Sw = cell(print_param.levels, 1); 
sol.traction = cell(print_param.levels, 1); 
sol.flux = cell(print_param.levels, 1); 

%% Flux controlled time loop
solver_param = struct(); % initialize structure to store solver parameters
solver_param.tol = 1E-8; % tolerance
solver_param.maxIter = 10; % maximum number of iterations
pControlled = false; % boolean to check bc control status
p = p_init; % current value of pressure
p_act = p; % actual value of pressure
u = u_init; % current value of displacement
p_top = min(p_init); % initial top pressure value

while (time_param.time < time_param.simTime) && (p_top > p_crit) ...
        && (pControlled == false)
      
    p_n = p; % current time level (n-index) 
    u_n = u; % current time level (n-index)
    time_param.time = time_param.time + time_param.tau; % cumulative time
    
    % Source terms
    sourceFlow = zeros(Nc, 1); % no sources for the flow
    sourceMech = modelEqsFlux.body(p_n); % sourceMech = body force
    
    % Calling Newton solver
    [p, p_m, u, iter] = solverUnsatBiot(G, p_n, u_n, modelEqsFlux, ...
        time_param, solver_param, sourceFlow, sourceMech);
    
    % Approximating top pressure
    fluxTemp = modelEqsFlux.Q(p, p_m);
    p_top = computeTopPressure(G, phys, p, fluxTemp, modelEqsFlux);
    
    % If it is flux controlled, update time step and store solution
    if (p_top > p_crit)
        
        % Calling time stepping routine
        [time_param.tau, print_param.print] = timeStepping(time_param, ...
            print_param, iter);
        
        % Actual pressure
        p_act = p;
        
        % Storing solutions
        if time_param.time == print_param.times(print_param.export)
            sol.time(print_param.export,1) = time_param.time;
            sol.u{print_param.export,1} = u;
            sol.p{print_param.export,1} = p_act;
            sol.pTop(print_param.export,1) = p_top;
            sol.Sw{print_param.export,1} = modelEqsFlux.S_w(p_act);
            sol.traction{print_param.export,1} = modelEqsFlux.T(u);
            sol.flux{print_param.export,1} = modelEqsFlux.Q(p_act, p_m);
            print_param.export = print_param.export + 1;
        end
        
    % If we reach the critical pressure, change to pressure controlled
    else
        fprintf('\n Changing from flux to pressure bc at %.1f [s]  \n\n', ...
            time_param.time - time_param.tau);
        pause(1);
        if print_param.print > 1 && ...
                print_param.times(print_param.print-1) == time_param.time
            print_param.print = print_param.print - 1;
        end
        time_param.time = time_param.time - time_param.tau; % back one step
        pControlled = true;
        time_param.tau = time_param.tau_min;
    end
end

%% Pressure controlled time loop
while (time_param.time < time_param.simTime) && (pControlled == true)
    
    p_n = p_act; % current time level (n-index) 
    u_n = u; % current time level (n-index)
    time_param.time = time_param.time + time_param.tau; % cumulative time
    
    % Source terms
    sourceFlow = zeros(Nc, 1); % no sources for the flow
    sourceMech = modelEqsPres.body(p_n); % sourceMech = body force
    
    % Calling Newton solver
    [p_act, p_m, u, iter] = solverUnsatBiot(G, p_n, u_n, modelEqsPres, ...
        time_param, solver_param, sourceFlow, sourceMech);
    
    % Approximating top pressure
    fluxTemp = modelEqsPres.Q(p_act, p_m);
    p_top = computeTopPressure(G, phys, p_act, fluxTemp, modelEqsPres);

    % Calling time stepping routine
    [time_param.tau, print_param.print] = timeStepping(time_param, ...
        print_param, iter);
    
    % Storing solutions
    if time_param.time == print_param.times(print_param.export)
        sol.time(print_param.export,1) = time_param.time;
        sol.u{print_param.export,1} = u;
        sol.p{print_param.export,1} = p_act;
        sol.pTop(print_param.export,1) = p_crit;
        sol.Sw{print_param.export,1} = modelEqsPres.S_w(p_act);
        sol.traction{print_param.export,1} = modelEqsPres.T(u);
        sol.flux{print_param.export,1} = modelEqsPres.Q(p, p_m);
        print_param.export = print_param.export + 1;
    end
    
end

% Displaying "sol" structure
fprintf('sol\n'); disp(sol);

%% Pressure plot
clf; minP = min(sol.p{end,1}); out = 1; 
for ii=1:length(sol.time)
    newplot; 
    plotCellData(G, sol.p{ii,1}); % plotting data/creating handle
    t_h=handle(text('Units','normalized', ... % Use [%] of the axis length
                    'Position',[.71,1.05,1], ... % Position of text   
                    'EdgeColor','k', ... % Textbox 
                    'FontSize', 13)); % Font size 
    t_h.String=sprintf('Time [hours]: %2.2f',sol.time(ii)/hour); % print time counter             
    axis equal off; view([-171 36]); cb = colorbar('South'); cb.FontSize=13; caxis([minP 0])
    title('p_w [Pa]','FontSize',15); set(gcf,'color','w');
    pause(0.01); out = out + 1;
    t_h.String = []; % delete the current text box to plot the next one
end

%% Saturation plot
clf; maxS = 1; minS = phys.flow.S_r; out = 1; 
for ii=1:length(sol.time)
    newplot;
    plotCellData(G, sol.Sw{ii,1}); % plotting data/creating handle
    t_h=handle(text('Units','normalized', ... % Use [%] of the axis length
                    'Position',[.71,1.05,1], ... % Position of text   
                    'EdgeColor','k', ... % Textbox 
                    'FontSize', 13)); % Font size 
    t_h.String=sprintf('Time [hours]: %2.2f',sol.time(ii)/hour); % print time counter             
    axis equal off; view([-171 36]); cb=colorbar('South'); cb.FontSize=13; caxis([minS maxS])
    title('S_w [ - ]','FontSize',15); set(gcf,'color','w');
    pause(0.01); out = out + 1; 
    t_h.String = []; % delete the current text box to plot the next one
end

%% Displacement plot
clf; u_mag = cell(length(sol.time), 1);
for ii=1:length(sol.time)
    u_mag{ii,1} = sqrt(sol.u{ii,1}(1:Nd:end).^2 ...
        + sol.u{ii,1}(2:Nd:end).^2 + sol.u{ii,1}(3:Nd:end).^2);
end
maxU = max(u_mag{end,1})/milli; minU = 0; out = 1; newplot;
for ii=1:length(sol.time)
    newplot;
    plotCellData(G, u_mag{ii,1}/milli); % plotting data/creating handle
    t_h=handle(text('Units','normalized',... % Use [%] of the axis length
                     'Position',[.75,1.05,1],... % Position of text   
                     'EdgeColor','k',... % Textbox 
                     'FontSize', 13)); % Font size 
    t_h.String=sprintf('Time [hours]: %2.2f',sol.time(ii)/hour); % print "time" counter             
    xlabel('x [m]'); ylabel('y [m]'); zlabel('Depth [m]');
    axis equal off; view([-171 36]); cb=colorbar('South'); cb.FontSize=13; caxis([minU maxU])
    title('|| u || [mm]','FontSize',15); box on; set(gcf,'color','w');
    pause(0.01); out = out + 1;
    t_h.String = []; % delete the current text box to plot the next one
end

%% Top pressure and fluxes plots
clf; newplot; 
fluxTopPlot = zeros(length(sol.time),1);
psiTopPlot = zeros(length(sol.time),1);
for ii=1:length(sol.time)
    fluxTopPlot(ii,1) = sum(abs(sol.flux{ii,1}(z_min)) / milli^3);
    psiTopPlot(ii,1) = (sol.pTop(ii,1)/phys.flow.gamma);
end
subplot(1,2,1);
plot(sol.time/hour,psiTopPlot,'r.-','LineWidth',2,'MarkerSize',15);
axis tight; grid on; box on;
xlabel('Time [hours]'); a = get(gca,'XTickLabel'); set(gca,'XTickLabel',a,'FontSize',14)
ylabel('\psi_w^{top} [m]  (x10^4)'); b = get(gca,'YTickLabel'); set(gca,'YTickLabel',b,'FontSize',14)
subplot(1,2,2);
plot(sol.time/hour,fluxTopPlot,'r.-','LineWidth',2,'MarkerSize',15);
axis tight; grid on; box on;
xlabel('Time [hours]'); a = get(gca,'XTickLabel'); set(gca,'XTickLabel',a,'FontSize',14)
ylabel('Q_w^{top}  [mm^3/s]'); b = get(gca,'YTickLabel'); set(gca,'YTickLabel',b,'FontSize',14)

%% Quiver plot of displacement field
clf; newplot;
ux = sol.u{end}(1:Nd:end); uxTop = ux(1:G.layerSize);
uy = sol.u{end}(2:Nd:end); uyTop = uy(1:G.layerSize);
plotCellData(G, u_mag{end,1}/milli, 1:G.layerSize); 
axis equal off; view([0 90]); colormap('gray');
cb=colorbar('EastOutside'); cb.FontSize=13; caxis([minU maxU])
title('|| u || [mm]','FontSize',15); box on; set(gcf,'color','w'); hold on;
quiver(xc(1:G.layerSize), yc(1:G.layerSize), uxTop, uyTop, 1,'LineWidth',1.5, 'Color', 'r');
pause(0.01); axis([0.05, 0.103, 0.05, 0.103]);