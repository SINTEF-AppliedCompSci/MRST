%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               fv-unsat                                                
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
clear; clc(); mrstModule add fvbiot upr

%% Creating grid

% Two-dimensional grid
r = 0.05;                         % [m] radii of the petri-dish
fd =@(p) sqrt(sum(p.^2, 2)) - r;  % circular domain function
min_x = -r;  max_x = r;           % [m] min and max values in x-axis
min_y = -r;  max_y = r;           % [m] min and max values in y-axis                  
h = (2*r)/25;                     % [m] step size
[p, t] = distmesh2d(fd, @huniform, h, [min_x, min_y; max_x, max_y], []);
p = p + r;                        % shifting triangulation points
G = triangleGrid(p, t);           % creating triangular grid
G = pebi(G);                      % creaing Voronoi diagram

% Extrude in the z-direction
Lx = max(p(:,1)); % max lenght x-axis
Ly = max(p(:,2)); % max lenght y-axis
Lz = 0.015; % [m] Thickness of the Petri-dish
nz = 5; % number of layers in z-axis
dz = Lz/nz; % thickness of each layer
thick = dz .* ones(nz,1);      % thickness vector
G = makeLayeredGrid(G, thick); % extrude grid
G = computeGeometry(G);        % compute geometry

%Plotting grid
newplot; plotGrid(G, 'FaceColor', [.8, .8, .8]); 
axis off equal; view([-9, 33]); 

%% Extracting grid information

Nc = G.cells.num;     % total number of cells
Nf = G.faces.num;     % total number of faces
Nd = G.griddim;       % grid dimension
V  = G.cells.volumes; % cell volumes
A  = G.faces.areas;   % face areas

xc = G.cells.centroids(:,1); % cell centers in x-direction
yc = G.cells.centroids(:,2); % cell centers in y-direction
zc = G.cells.centroids(:,3); % cell centers in z-direction
        
xf = G.faces.centroids(:,1); % face centers in x-direction
yf = G.faces.centroids(:,2); % face centers in y-direction
zf = G.faces.centroids(:,3); % face centers in z-direction

nux = G.faces.normals(:,1); % face normals in x-direction
nuy = G.faces.normals(:,2); % face normals in y-direction
nuz = G.faces.normals(:,3); % face normals in z-direction

zetac = Lz - zc; % cell centers of elev. head
zetaf = Lz - zf; % face centers of elev. head

fNei = G.faces.neighbors; % extracting faces neighbors
int_f = find(all(fNei ~= 0,2)); % internal faces
ext_f = find(((fNei(:,1) == 0) + (fNei(:,2) == 0)) == 1); % external faces
int_fNei = fNei(all(fNei ~= 0,2),:); % internal faces neighbors

z_min = find(zf == 0); % idx of top faces
z_max = find(zf > 0.9999*Lz & zf < 1.0001*Lz ); % idx of bottom faces
sides = ext_f(~ismember(ext_f, [z_min; z_max])); % idx of sides faces

%% Declaring physical parameters

% Mechanics parameters (Kaolinite (BS) Table 2 -> Mondol et al, 2008)
lambda_s = 1.229E11 .* ones(Nc, 1);  % [Pa] first Lame parameter
mu_s = 4.7794E10 .* ones(Nc, 1);     % [Pa] second Lame parameter
C_s = 5.618E-11;                     % [1/Pa] solid compressibility
rho_s = 1769;                        % [kg/m^3] solid density
stiff = shear_normal_stress(Nc, Nd, mu_s, lambda_s, 0.*mu_s); % stiff mat

% Flow parameters
soil = getHydraulicProperties('clay');
K_sat = soil.K_s;            % [m/s] saturated hydraulic conduct.
rho_w = 1014;                % [kg/m^3] water density
g = 9.8006;                  % [m/s^2] gravity acceleration
gamma = rho_w * g;           % [Pa/m] specific gravity
mu_w = 1 * centi * poise;    % [Pa*s] water dynamic viscosity
k = K_sat*mu_w/gamma;        % [m^2] intrinsic permeability
rock.perm = k * ones(Nc, 1); % permeability structure
n = soil.theta_s;            % [-] porosity
C_w = 4.5455E-10;            % [1/Pa] water compressibility
alphaVan = soil.alpha;       % [1/m] Equation parameter
aVan = alphaVan / gamma;     % [Pa/m] Artificial parameter
nVan = soil.n;               % [-] Equation parameter
mVan = 1 - ( 1 / nVan);      % [-] Equation parameter
theta_r = soil.theta_r;      % [-] Residual water content
S_r = theta_r / n;           % [-] Residual saturation
tempAmbient = 20;            % [C] Ambient temperature
relHumidity = 50;            % [%] Ambient relative humidity
p_crit = computeCriticalPressure(tempAmbient, relHumidity, rho_w);

% Coupling parameters
C_m = 2.17E-10;               % [1/Pa] porous medium compressibility
alpha_biot = 1 - C_s/C_m;     % [-] Biot's coefficient

%% Calling water rentention model

[S_w, krw, C_S] = vanGenuchtenMualemSw(aVan, S_r, nVan, mVan);

%% Boundary conditions

% MECHANICS BOUNDARY CONDITIONS
bcMech = addBC([], sides, 'pressure', 0);     % u=0 at the sides
bcMech = addBC(bcMech, z_max, 'pressure', 0); % u=0 at the bottom

% Boundary values are assigned to a vector of size Nd * Nf. Face
% number one will have its x-condition assigned to bc_vals(1), y-condition
% at bc_vals(2), and z-condition at bc_vals(3). Next face number 2 etc.
% Face i has traction (or displacement) index 3*i-2, 3*i-1, 3*i (3D).
% Mixed boundary conditions, i.e., roller boundary conditions are not
% implemented

bcMechVals = zeros(Nd*Nf, 1);

% Setting lateral faces as zero displacement 
bcMechVals(Nd*sides-2) = 0; % ux = 0 at the sides
bcMechVals(Nd*sides-1) = 0; % uy = 0 at the sides
bcMechVals(Nd*sides) = 0;   % uz = 0 at the sides

% Setting bottom layer as zero displacement
bcMechVals(Nd*z_max-2) = 0; % ux = 0 at the bottom
bcMechVals(Nd*z_max-1) = 0; % uy = 0 at the bottom
bcMechVals(Nd*z_max) = 0;   % uz = 0 at the bottom

% FLOW BOUNDARY CONDITONS

% FLUX CONTROLLED BOUNDARY CONDITIONS
topArea = sum(A(z_min));  % [m^2] surface top area
vMax = 6E-07;             % [m/s] maximum evaporation velocity
EMax = vMax * topArea;    % [m^3/s] maximum evaporation flux

% Now we create a vector containing the indidual fluxes for each
% top face. Note that we mulitply by the viscosity, since the 
% treatment of the boundary conditions for the mpfa routine assumes
% unit viscosity.
Qtop_f = vMax .* mu_w .* A(z_min) ; % [m^3/s] Top boundary

% Creating the boundary structure for flux controlled BC
bcFlow_f            = addBC([], z_min, 'flux', Qtop_f);     
bcFlowVals_f        = zeros(Nf, 1); 
bcFlowVals_f(z_min) = Qtop_f;                        

% PRESSURE CONTROLLED BOUNDARY CONDITIONS

% Creating the boundary structure for pressure controlled BC
bcFlow_p            = addBC([], z_min, 'pressure', p_crit);  
bcFlowVals_p        = zeros(Nf, 1);                          
bcFlowVals_p(z_min) = p_crit + gamma .* zetaf(z_min);                             
% The second term of the above equation represents the gravity 
% contribution, and must be included for Dirichlet boundary conditions

%% Calling MPSA/MPFA routines and creating discrete operators

% Discretize mechanics problem
mpsa_discr = mpsa(G, stiff, [], 'invertBlocks', 'matlab', 'bc', bcMech);

% Discretize flow problem for flux-controlled boundary conditions
mpfa_discr_flux = mpfa(G, rock, [], 'invertBlocks', 'matlab', 'bc', bcFlow_f);

% Discretize flow problem for pressure-controlled boundary conditions
mpfa_discr_pres = mpfa(G, rock, [], 'invertBlocks', 'matlab', 'bc', bcFlow_p);

% Mechanics discrete operators
S           = @(x) mpsa_discr.stress * x; % stres
boundS      = @(x) mpsa_discr.boundStress * x; % stress boundary
gradP       = @(x) mpsa_discr.gradP * x; % gradient of pressure
divU        = @(x) mpsa_discr.divD * x; % divergence of displacement
compat      = @(x) mpsa_discr.stabDelta * x; % stability operator
divS        = @(x) mpsa_discr.div * x; % divergence of stress

% Flow discrete operators
F_f         = @(x) mpfa_discr_flux.F * x; % flux (flux-control)
F_p         = @(x) mpfa_discr_pres.F * x; % flux (p-control)
boundF_f    = @(x) mpfa_discr_flux.boundFlux * x; % bcflux (flux-control)
boundF_p    = @(x) mpfa_discr_pres.boundFlux * x; % bocflux (p-control)
divF        = @(x) mpfa_discr_flux.div * x; % divergence of flux

%% Creating AD variables
u_init = zeros(Nd*Nc,1);            % initial u
p_init = -0.1 * kilo *  ones(Nc,1); % initial p
u_ad = initVariablesADI(u_init);    % initializing u_ad
p_ad = initVariablesADI(p_init);    % initializing p_ad

%% Time parameters and printing parameters
simTime = 2*hour;       % [s] final simulation time
tau_init = 25;          % [s] initial time step
tau_min = 25;           % [s] minimum time step
tau_max = 1000;         % [s] maximum time step
tau = tau_init;         % [s] initializing time step
timeCum = 0;            % [s] initializing cumulative time

printLevels = 50;       % number of printing levels
printTimes = ((simTime/printLevels):(simTime/printLevels):simTime)';
pp = 1;                 % initializing printing counter
ee = 1;                 % intializing exporting counter

%% Discrete equations

% Function that maps "c" to "Nd*c"
sca2vec = @(x) repmat(x, [Nd, 1]); 

% MECHANICS DISCRETE EQUATIONS

% Body forces
g_body = zeros(Nd*Nc,1); g_body(Nd:Nd:end) = g; % assigning g to z-cells
body = @(p_n) ((1-n).*rho_s + n.*sca2vec(S_w(p_n)).*rho_w) .* g_body;
% Traction
T = @(u) S(u) + boundS(bcMechVals);
% Momentum equation (Mechanics contribution)
uEq1 = @(u) divS(T(u));
% Momentum equation (Flow contribution + source)
uEq2 = @(p, p_n, sourceMech) -alpha_biot .* gradP(S_w(p_n) .* p) ...
    + sca2vec(V) .* sourceMech;

% FLOW DISCRETE EQUATIONS

% Flux-controlled bc discrete equations

% Upstream weighting of krw
krwUp_f = @(p_m) upstreamWeightingMPFA(G, bcFlow_f, bcFlowVals_f, ...
    mpfa_discr_flux, krw, gamma, p_m, 'pressure', 'on');
% Compressibility-like terms
xi  = @(p_n) (alpha_biot-n).*C_s.*S_w(p_n).^2 + n.*C_w.*S_w(p_n);
chi = @(p_n) (alpha_biot-n).*C_s.*S_w(p_n).*p_n + n;  
% Darcy Flux
Q_f = @(p, p_m) (krwUp_f(p_m) ./ mu_w) .* (F_f(p + gamma.*zetac) ...
    + boundF_f(bcFlowVals_f)); 
% Mass conservation (Mechanics contribution)
pEq1 = @(p_n, u, u_n) alpha_biot .* S_w(p_n) .* divU(u-u_n);
% Mass conservation (Flow contribution)
pEq2_f = @(p, p_n, p_m, tau, sourceFlow)  ...
       alpha_biot.^2 .* S_w(p_n) .* compat(S_w(p_n).*p_n) ...
       + V .* xi(p_n) .* (p - p_n) ...
       + V .* chi(p_n) .* (S_w(p_m) + C_S(p_m) .* (p - p_m) - S_w(p_n)) ...
       + tau .* divF(Q_f(p, p_m)) - V .* sourceFlow;
     
% Pressure-controlled bc discrete equations

% Upstream Weighting of krw
krwUp_p = @(p_m) upstreamWeightingMPFA(G, bcFlow_p, bcFlowVals_p, ...
    mpfa_discr_pres, krw, gamma, p_m, 'pressure', 'on');
% Darcy Flux
Q_p = @(p, p_m) (krwUp_p(p_m) ./ mu_w) .* (F_p(p + gamma.*zetac) ...
    + boundF_p(bcFlowVals_p));
% Mass conservation (Flow contribution)
pEq2_p = @(p, p_n, p_m, tau, sourceFlow)  ...
       alpha_biot.^2 .* S_w(p_n) .* compat(S_w(p_n).*p_n)...
       + V .* xi(p_n) .* (p - p_n) ...
       + V .* chi(p_n) .* (S_w(p_m) + C_S(p_m) .* (p - p_m) - S_w(p_n)) ...
       + tau .* divF(Q_p(p, p_m)) - V .* sourceFlow;   
 
%% Solution structure
sol.time      = zeros(printLevels, 1); % time
sol.u         = cell(printLevels, 1);  % displacement
sol.p         = cell(printLevels, 1);  % pressure
sol.pTop      = zeros(printLevels, 1); % average top pressure
sol.Sw        = cell(printLevels, 1);  % saturation
sol.traction  = cell(printLevels, 1);  % traction
sol.flux      = cell(printLevels, 1);  % flux

%% Flux controlled time loop

pControlled = false;             % boolean to check bc control status
p_act = p_ad;                    % actual value of pressure
p_top = min(p_init);             % intial top pressure value
tol = 1E-8;                      % tolerance
maxIter = 10;                    % maximum number of iterations
while (timeCum < simTime) && (p_top > p_crit) && (pControlled == false)
      
    p_n = p_ad.val; % current time level (n-index) 
    u_n = u_ad.val; % current time level (n-index)
    timeCum = timeCum + tau; % cumulative time
    
    % Source terms
    sourceFlow = zeros(Nc, 1);  % no sources for the flow
    sourceMech = body(p_n);     % sourceMech = body force
    
    % Calling Newton solver
    [p_ad, p_m, u_ad, iter] = solverUnsatBiot(G, p_ad, p_n, ...
        u_ad, u_n, pEq1, pEq2_f, uEq1, uEq2, tau, sourceFlow, sourceMech, ...
        timeCum, tol, maxIter);
    
    % Approximating top pressure
    fluxTemp = Q_f(p_ad.val, p_m); krwTemp = krwUp_f(p_m);
    p_top = computeTopPressure(G, p_ad, fluxTemp, gamma, mu_w, k, krw);
    
    % If it is flux controlled, update time step and store solution
    if (p_top > p_crit)
        
        % Calling time stepping routine
        [tau, pp] = timeStepping(tau, tau_min, tau_max, simTime, timeCum,...
        iter, printTimes, pp);
        
        % Actual pressure
        p_act = p_ad;
        
        % Storing solutions
        if timeCum == printTimes(ee)
            sol.time(ee,1)      = timeCum;
            sol.u{ee,1}         = u_ad.val;
            sol.p{ee,1}         = p_act.val;
            sol.pTop(ee,1)      = p_top;
            sol.Sw{ee,1}        = S_w(p_act.val);
            sol.traction{ee,1}  = T(u_ad.val);
            sol.flux{ee,1}      = Q_f(p_act.val, p_m);
            ee = ee + 1;
        end
        
    % If we reach the critical pressure, change to pressure controlled
    else
        fprintf('\n Changing from flux to pressure bc at %.1f [s]  \n\n', ...
            timeCum-tau);
        pause(1);
        if pp > 1 && printTimes(pp-1) == timeCum
            pp = pp - 1;
        end
        timeCum = timeCum-tau; % we go back one time step
        pControlled = true;
        tau = tau_min;
    end
end

%% Pressure controlled time loop
while (timeCum < simTime) && (pControlled == true)
    
    p_n = p_act.val; % current time level (n-index) 
    u_n = u_ad.val;  % current time level (n-index)
    timeCum = timeCum + tau;     % cumulative time
    
    % Source terms
    sourceFlow = zeros(Nc, 1);  % no sources for the flow
    sourceMech = body(p_n);     % sourceMech = body force
    
    % Calling Newton solver
    [p_act, p_m, u_ad, iter] = solverUnsatBiot(G, p_act, p_n, ...
        u_ad, u_n, pEq1, pEq2_p, uEq1, uEq2, tau, sourceFlow, sourceMech, ...
        timeCum, tol, maxIter);
    
    fluxTemp = Q_p(p_act.val, p_m); 
    krwTemp = krwUp_p(p_m);
    p_top = computeTopPressure(G, p_act, fluxTemp, gamma, mu_w, k, krw);
       
    % Calling time stepping routine
    [tau, pp] = timeStepping(tau, tau_min, tau_max, simTime, timeCum,...
        iter, printTimes, pp);
    
    % Storing solutions
    if timeCum == printTimes(ee)
        sol.time(ee,1)      = timeCum;
        sol.u{ee,1}         = u_ad.val;
        sol.p{ee,1}         = p_act.val;
        sol.pTop(ee,1)      = p_crit;
        sol.Sw{ee,1}        = S_w(p_act.val);
        sol.traction{ee,1}  = T(u_ad.val);
        sol.flux{ee,1}      = Q_p(p_act.val, p_m);
        ee = ee + 1;
    end
    
end

% Storing final solutions in "sol" structure
fprintf('sol\n'); disp(sol);

%% Pressure plot
clf; maxP = max(sol.p{1,1}); minP = min(sol.p{end,1}); out = 1; 
for ii=1:length(sol.time)
    newplot; 
    plotCellData(G, sol.p{ii,1});               % plotting data/creating handle
    t_h=handle(text('Units','normalized',...    % Use [%] of the axis length
                    'Position',[.71,1.05,1],... % Position of text   
                    'EdgeColor','k',...         % Textbox 
                    'FontSize', 13));           % Font size 
    t_h.String=sprintf('Time [hours]: %2.2f',sol.time(ii)/hour); % print time counter             
    axis equal off; view([-171 36]); cb = colorbar('South'); cb.FontSize=13; caxis([minP 0])
    title('p_w [Pa]','FontSize',15); set(gcf,'color','w');
    pause(0.01); out = out + 1;
    t_h.String = []; % we delete the current text box to plot the next one
end

%% Saturation plot
clf; maxS = 1; minS = S_r; out= 1; 
for ii=1:length(sol.time)
    newplot;
    plotCellData(G, sol.Sw{ii,1});              % plotting data/creating handle
    t_h=handle(text('Units','normalized',...    % Use [%] of the axis length
                    'Position',[.71,1.05,1],... % Position of text   
                    'EdgeColor','k',...         % Textbox 
                    'FontSize', 13));           % Font size 
    t_h.String=sprintf('Time [hours]: %2.2f',sol.time(ii)/hour); % print time counter             
    axis equal off; view([-171 36]); cb=colorbar('South'); cb.FontSize=13; caxis([minS maxS])
    title('S_w [ - ]','FontSize',15); set(gcf,'color','w');
    pause(0.01); out = out + 1; 
    t_h.String = [];  % we delete the current text box to plot the next one
end

%% Displacement plot
clf; u_mag = cell(length(sol.time), 1);
for ii=1:length(sol.time)
    u_mag{ii,1} = sqrt(sol.u{ii,1}(1:Nd:end).^2 ...
        + sol.u{ii,1}(2:Nd:end).^2 + sol.u{ii,1}(3:Nd:end).^2);
end
maxU = max(u_mag{end,1})/milli; minU = 0; out= 1; newplot;
for ii=1:length(sol.time)
    newplot;
    plotCellData(G, u_mag{ii,1}/milli);          % plotting data/creating handle
    t_h=handle(text('Units','normalized',...     % Use [%] of the axis length
                     'Position',[.75,1.05,1],... % Position of text   
                     'EdgeColor','k',...         % Textbox 
                     'FontSize', 13));            % Font size 
    t_h.String=sprintf('Time [hours]: %2.2f',sol.time(ii)/hour); % print "time" counter             
    xlabel('x [m]'); ylabel('y [m]'); zlabel('Depth [m]');
    axis equal off; view([-171 36]); cb=colorbar('South'); cb.FontSize=13; caxis([minU maxU])
    title('|| u || [mm]','FontSize',15); box on; set(gcf,'color','w');
    pause(0.01);    out = out + 1;
    t_h.String = [];   % we delete the current text box to plot the next one
end

%% Top pressure and fluxes plots
clf; newplot;
fluxTopPlot = zeros(length(sol.time),1);
psiTopPlot = zeros(length(sol.time),1);
for ii=1:length(sol.time)
    fluxTopPlot(ii,1) = sum(abs(sol.flux{ii,1}(z_min)) / milli^3);
    psiTopPlot(ii,1) = (sol.pTop(ii,1)/gamma);
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

%% Quiver plot
clf; newplot;
ux = sol.u{end}(1:Nd:end); uxTop = ux(1:G.layerSize);
uy = sol.u{end}(2:Nd:end); uyTop = uy(1:G.layerSize);
plotCellData(G, u_mag{end,1}/milli, 1:G.layerSize); % plotting data/creating handle
axis equal off; view([0 90]); colormap('gray');
cb=colorbar('EastOutside'); cb.FontSize=13; caxis([minU maxU])
title('|| u || [mm]','FontSize',15); box on; set(gcf,'color','w'); hold on;
quiver(xc(1:G.layerSize), yc(1:G.layerSize), uxTop, uyTop, 1,'LineWidth',1.5, 'Color', 'r');
pause(0.01); axis([0.05, 0.103, 0.05, 0.103]);