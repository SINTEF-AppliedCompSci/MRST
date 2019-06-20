%% Three-dimensional water infiltration in a heterogeneous domain 
% 
% This code is an implementation of the mixed-based form of the Richards' 
% Equation in a three-dimensional domain using Automatic Differentiation 
% from MRST (see <http://www.sintef.no/projectweb/mrst/>). 
% The grid is set to be cartesian and structured. For the spatial
% discretization we use cell-centered finite volume method with multi-point
% flux approximation from <https://github.com/keileg/fvbiot>. For the time 
% derivative we use the modified picard iteration proposed in: 
% http://onlinelibrary.wiley.com/doi/10.1029/92WR01488/abstract.
% We decided to use the hydraulic head as primary variable instead of 
% the classical approach (the pressure head) in order to avoid any
% inconsistency in the gravity contribution.
% The intrisic permeability at the faces are approximated using harmonic 
% average while the relative permeabilities are upstream weighted. The 
% problem of interest is the vertical infiltration of water 
% (from top to bottom). Moreover, we included low permeability formations
% near the top and bottom boundaries.
% For the boundary conditions we use constant $\psi$ at top and bottom 
% of the domain and no flux elsewhere.
%
% *Equations:*
%
% Mass Conservation
%
% $$ \frac{\partial \theta}{\partial t} + 
%     \nabla \cdot \vec{q} = 0 $$
%                 
% Multiphase Darcy's Law
%
% $$ \vec{q} = -k \frac{\rho g}{\mu} k_{rw} 
%             \nabla h $$
%
% *Discrete equations:*
%
% Mass Conservation
%
% $$ \frac{V_{\{c\}}}{\Delta t} \left( 
%    \theta_{\{c\}}^{n+1,m}  + 
%     C_{\{c\}}^{n+1,m} 
%    \left( h_{\{c\}}^{n+1,m+1} - h_{\{c\}}^{n+1,m} \right) 
%    - \theta_{\{c\}}^n \right) 
%    + \left[ \mathbf{div} \left( \vec{q}_{\{f\}} \right) \right]_{\{c\}} = 0 $$
%
% Multiphase Darcy's Law
%
% $$ \vec{q}_{\{f\}} = \frac{\rho g}{\mu} 
%    \left[ \mathbf{upstr} \left( {k_{rw}}_{\{c\}}^{n+1,m} \right) \right]_{\{f\}}  
%    \left( \left[ \mathbf{F} \left( h_{\{c\}}^{n+1,m+1}\right)\right]_{\{f\}} 
%    + \left[ \mathbf{boundF} \left(b_{\{f\}} \right) \right]_{\{f\}} \right)$$
%
% *Operators, variables and parameters*
% 
% $\vec{\chi}$ : vector quantity
% 
% $\nabla \chi$ : gradient operator
% 
% $\nabla \cdot \chi$ : divergence operator
%
% $\chi_{\{c\}}$ : cell evaluated quantity
%
% $\chi_{\{f\}}$ : face evaluated quantity
%
% $\chi^{n}$ : quantity evaluated at the current time step
%
% $\chi^{n+1}$ : quantity evaluated at the next time step
%
% $\chi^{m}$ : quantity evaluated at the current iteration level
%
% $\chi^{m+1}$ : quantity evaluated at the next iteration level
%
% $\mathbf{div} \left( \chi \right)$ : mpfa divergence operator (maps faces
% to centers).
%
% $\mathbf{upstr} \left( \chi \right)$ : upstream weighted quantity (maps centers to
% faces by determining the upstream flux).
% 
% $\mathbf{F} \left( \chi \right)$ : mpfa discrete gradient operator (maps centers to
% faces).
%
% $\mathbf{boundF} \left( \chi \right)$ : operator that deals with
% discretization of boundary values (maps faces to faces).
%
% $h := \psi + \zeta$ : hydraulic head $[L]$ (see <Pinder, G. F., & Celia, M.
% A. (2006). Subsurface hydrology>)
%
% $\psi$ : pressure head $[L]$
%
% $\zeta$ : elevation head $[L]$
%
% $\theta = \theta(\psi)$ : water content $[-]$
%
% $v$ : water Darcy's velocity $[L/T]$
%
% $k$ : intrisic permeability $[L^2]$
%
% $\rho$ : water density $[M/L^3]$
%
% $g$ : acceleration of gravity $[L/T^2]$
%
% $\mu$ : water dynamic viscosity $[M/LT]$
%
% $k_{rw} = k_{rw}(\psi)$ : water relative permeability $[-]$
%
% $V$ : volume of REV $[L^3]$
%
% $\Delta t$ : time step length $[T]$
%
% $C = C(\psi)$ : specific moisture capacity $[1/L]$
% 
% $q$ : water Darcy's flux $[L^3/T]$
%
% $b$ : boundary values vector ($[L]$ if it's hydraulic head, $[L^3/T]$ if it's
% flux).
%
% *van Genuchten parameters*
%
% $$ \theta (\psi) = \frac{\theta_s - \theta_r}{\left[1 + (\alpha 
% |\psi |)^n \right]^m} $$,
%
% $$ C(\psi) = \frac{d\theta}{d\psi} = \frac{mn\psi (-(\theta_s-\theta_r))
%        \alpha^n |\psi|^{n-2}}{(\alpha^n |\psi|^n +1)^{m+1}} $$,
%
% $$ k_r(\psi) = \frac{ \big\{ 1 - (\alpha |\psi|^{n-1}) 
% [1 +  (\alpha |\psi|^n)]^{-m}  \big\}^2}{[1+(\alpha |\psi|^n)]^{m/2}} $$,
%
% where $\theta_s$ is the soil moisture content, $\theta_r$ is
% the residual moisture content and $\alpha$, $n$ and $m$ are 
% the van Genuchten's parameters. 
%
% *Physical Parameters* (Taken from:
% http://onlinelibrary.wiley.com/doi/10.1029/92WR01488/abstract)
% 
% * $\alpha = 0.0335 \; [1/cm]$ 
%
% * $\theta_s = 0.368$ 
%
% * $\theta_r = 0.102$
%
% * $n = 2$
%
% * $m = 0.5$
%
% * $K_{sat} = 0.00922 \, [cm \cdot s^{-1}]$  (sat. hydraulic conductivity)
% 
% * $\rho = 1 \; [g\cdot cm^{-3}]$
%
% * $\mu = 1 \; [cP]$ 
%
% * $g = 981 \; [cm\cdot s^{-2}]$
%
% *Boundary Conditions*
% 
% * Top boundary: $\psi = - 75 \; [cm]$ 
%
% * Bottom boundary: $\psi = -1000 \; [cm]$
% 
% *Intial Conditions*
%
% * $\psi = -1000 \; [cm]$
%

%% Clearing workspace and cleaning console
clear; clc();

%% Importing modules
mrstModule add re-mpfa fvbiot

%% Setting up the Grid
nx = 10;              % Cells in x-direction
ny = 10;              % Cells in y-direction
nz = 40;              % Cells in z-direction
Lx = 100;             % Lenght in x-direction
Ly = 100;             % Length in y-direction
Lz = 100;             % Length in z-direction
G = computeGeometry(cartGrid([nx,ny,nz],[Lx,Ly,Lz])); % computing geometry
V = G.cells.volumes;  % Cell volumes

% Plotting Grid
figure(1); plotGrid(G);
xlabel('x-axis [cm]'); ylabel('y-axis [cm]'); zlabel('Depth [cm]');
axis tight; pbaspect([1,1,4]); view([-51,26]); 

%% Fluid Properties
rho = 1;                    % [g/cm^3] water density
mu = 0.01;                  % [g/cm.s] water viscosity
g = 980.6650;               % [cm/s^2] gravity acceleration

%% Rock Properties
K_sat = 0.00922;            % [cm/s] saturated hydraulic conductivity
k = (K_sat*mu)/(rho*g);     % [cm^2] intrinsic permeability
kz = 1E-20 * k;             % [cm^2] low permeability formation
rock.perm = repmat([k, k, k], [G.cells.num, 1]);  % creating perm structure

% Let's add a region of low permeability near the top boundary
rock.perm((12*nx*ny)+1:(12*nx*ny)+(nx*ny),3) = kz; 
rock.perm((12*nx*ny)+45:(12*nx*ny)+46,3) = k;
rock.perm((12*nx*ny)+55:(12*nx*ny)+56,3) = k;

% Let's add a region of low permeability near the bottom boundary
rock.perm((28*nx*ny)+1:(28*nx*ny)+(nx*ny),3) = kz; 
rock.perm((28*nx*ny)+45:(28*nx*ny)+46,3) = k;
rock.perm((28*nx*ny)+55:(28*nx*ny)+56,3) = k;

%% Plotting Grid and Low Permeability regions
subplot(1,2,1)
   plotGrid(G);
   xlabel('x-axis [cm]'); ylabel('y-axis [cm]'); zlabel('Depth [cm]');
   view([-40 12]); set(gcf,'color','w'); 
   title('Computational Grid');
subplot(1,2,2)
   plotCellData(G,rock.perm(:,3),rock.perm(:,3) < 1E-20);
   plotGrid(G, 'FaceAlpha', 0, 'EdgeAlpha', .1)
   xlabel('x-axis [cm]'), ylabel('y-axis [cm]'), zlabel('Depth [cm]');
   view([-40 12]); shading faceted; camproj perspective; axis tight;
   set(gca, 'ZDir', 'reverse'); box on;
   title('Low Permeability Regions');
      
%% Van Genuchten Parameters
alpha = 0.0335;       % [1/cm] Equation parameter
nVan = 2;             % [-]    Equation parameter
mVan = 1-(1/nVan);    % [-]    Equation parameter
theta_s = 0.368;      % [-]    Saturation soil moisture
theta_r = 0.102;      % [-]    Residual soil moisture

%% Boundary and Initial Conditions

% Extracting Grid information
zCentr = G.cells.centroids(:,3);   % centroids of cells in z-direction
zFaces = G.faces.centroids(:,3);   % centroids of faces in z-direction
zetaCentr = Lz - zCentr;           % centroids of cells of elev. head
zetaFaces = Lz - zFaces;           % centroids of faces of elev. head

% Determining boundary indices
x_min = find(G.faces.centroids(:,1) == 0);             % west bound indices
x_max = find(G.faces.centroids(:,1) > 0.9999*Lx & ...   
             G.faces.centroids(:,1) < 1.0001*Lx );     % east bound indices

y_min = find(G.faces.centroids(:,2) == 0);             % south bound indices
y_max = find(G.faces.centroids(:,2) > 0.9999*Ly & ...   
             G.faces.centroids(:,2) < 1.0001*Ly );     % north bound indices

z_min = find(G.faces.centroids(:,3) == 0);             % top bound indices
z_max = find( G.faces.centroids(:,3) > 0.9999*Lz & ...   
              G.faces.centroids(:,3) < 1.0001*Lz );    % bottom bound indices

% Boundary values          
fluxW = 0;                                  % [cm^3/s] West boundary 
fluxE = 0;                                  % [cm^3/s] East boundary   
fluxS = 0;                                  % [cm^3/s] South boundary
fluxN = 0;                                  % [cm^3/s] North boundary
psiT = -75;                                 % [cm]     Top boundary
psiB = -1000;                               % [cm]     Bottom boundary

hT = psiT + zetaFaces(z_min);               % [cm] Top hydraulic head
hB = psiB + zetaFaces(z_max);               % [cm] Bottom hydrualic head

% Creating the boundary structure
bc = addBC([], x_min, 'flux', fluxW);       % setting West boundary
bc = addBC(bc, x_max, 'flux', fluxE);       % setting East boundary    
bc = addBC(bc, y_min, 'flux', fluxS);       % setting South boundary
bc = addBC(bc, y_max, 'flux', fluxN);       % setting North boundary
bc = addBC(bc, z_min, 'pressure', hT);      % setting Top boundary
bc = addBC(bc, z_max, 'pressure', hB);      % setting Bottom boundary 

% Setting the boundary values vector
bc_val = zeros(G.faces.num, 1);             % initializing
bc_val(x_min) = fluxW;                      % assigning West boundary
bc_val(x_max) = fluxE;                      % assigning East boundary
bc_val(y_min) = fluxS;                      % assigning South boundary
bc_val(y_max) = fluxN;                      % assigning North boundary
bc_val(z_min) = hT;                         % assigning Top boundary
bc_val(z_max) = hB;                         % assigning Bottom boundary

% Initial Condition
h_init = psiB + zetaCentr;                  % [cm] Initial condition

%% Calling MPFA routine
mpfa_discr = mpfa(G,rock,[],'bc',bc,'invertBlocks','matlab');

%% Creating discrete mpfa-operators      
F = @(x) mpfa_discr.F * x;                  % flux discretization
boundF = @(x) mpfa_discr.boundFlux * x;     % boundary discretization
divf = @(x) mpfa_discr.div * x;             % divergence

%% Creating AD variable
h_ad = initVariablesADI(h_init);

%% Water retention curves

% Boolean function to determine if we are in the unsat or sat zone
isUnsat = @(x) x < 0;

% Water content
theta= @(x)     isUnsat(x) .* ((theta_s - theta_r) ./ ...
                (1 + (alpha .* abs(x)).^nVan).^mVan + theta_r ) + ...
                ~isUnsat(x) .* theta_s;
    
% Specific Moisture Capacity
cVan= @(x)      isUnsat(x) .* ((mVan .* nVan .* x .* (theta_r-theta_s) .* ...
                alpha.^nVan .* abs(x).^(nVan-2) ) ./ ...
                (alpha^nVan .* abs(x).^nVan + 1).^(mVan+1)) + ...
                ~isUnsat(x) .* 0;

% Relative permeability
krw= @(x)       isUnsat(x) .* ((1 - (alpha .* abs(x)).^(nVan-1) .* ...
                (1 + (alpha .* abs(x)).^nVan).^(-mVan)).^2 ./ ...
                (1 + (alpha .* abs(x)).^nVan).^(mVan./2)) + ...
                ~isUnsat(x) .* 1;

%% Time parameters
simTime = 96*3600;      % [s]    final simulation time
dt_init = 0.01;         % [s]    initial time step
dt_min = 0.01;          % [s]    minimum time step 
dt_max = 10000;         % [s]    maximum time step
lowerOptIterRange = 3;  % [iter] lower optimal iteration range
upperOptIterRange = 7;  % [iter] upper optimal iteration range
lowerMultFactor = 1.3;  % [-]    lower multiplication factor
upperMultFactor = 0.7;  % [-]    upper multiplication factor
dt = dt_init;           % [s]    initializing time step
timeCum = 0;            % [s]    initializing cumulative time
currentTime = 0;        % [s]    current time

%% Printing parameters
printLevels = 20;       % number of printing levels
printTimes = ((simTime/printLevels):(simTime/printLevels):simTime)';   % printing times
printCounter = 1;       % initializing printing counter
exportCounter = 1;      % intializing exporting counter

%% Discrete equations

% Upstream Weighting
krwUp = @(h_m0) upstreamMPFA_hydHead(krw,F,boundF,G,bc,bc_val,h_m0,...
                                     zetaCentr,zetaFaces);                  
                                  
% Darcy Flux
q = @(h,h_m0)  (rho.*g./mu) .* krwUp(h_m0) .* (F(h) + boundF(bc_val));

% Mass Conservation                         
hEq = @(h,h_n0,h_m0,dt) (V./dt) .* (...
                                theta(h_m0-zetaCentr) + ...
                                cVan(h_m0-zetaCentr) .* (h - h_m0) - ...
                                theta(h_n0-zetaCentr) ...
                                ) + divf(q(h,h_m0));
                            
%% Creating solution structure

sol = struct('time',[],'h',[],'psi',[],'zeta',[],'theta',[],'flux',[],...
             'iter',[],'residual',[],'cpuTime',[]);

% We need to initialize doubles and cells to store the values in sol   
t = zeros(printLevels,1);            % Cummulative time
h_cell = cell(printLevels,1);        % Hydraulic Heads
psi_cell = cell(printLevels,1);      % Pressure Heads
zeta_cell = cell(printLevels,1);     % Elevation Heads
theta_cell = cell(printLevels,1);    % Water content
flux_cell = cell(printLevels,1);     % Fluxes
residual = cell(printLevels,1);      % Residuals
iterations = zeros(printLevels,1);   % Iterations
cpuTime = zeros(printLevels,1);      % CPU time                                       
                            
%% Solving via Newton

while timeCum < simTime
    
    % Newton parameters
    maxTolPresHead = 1; % [cm] maximum absolute tolerance of pressure head
    maxIter = 10;       % maximum number iterations

    % Calling Newton solver
    [h_ad,h_m0,iter,timeCum,tf] = newtonSolver(h_ad,hEq,timeCum,...
                                               dt,maxTolPresHead,maxIter);
    
    % Calling Time stepping routine
    [dt,printCounter] = timeStepping(dt,dt_min,dt_max,simTime,...
                           timeCum,iter,printTimes,printCounter,...
                           lowerOptIterRange,upperOptIterRange,...
                           lowerMultFactor,upperMultFactor);
    
    % Storing solutions at each printing time
    if timeCum == printTimes(exportCounter)
        h_cell{exportCounter,1} = h_ad.val;
        psi_cell{exportCounter,1} = h_ad.val - zetaCentr;
        zeta_cell{exportCounter,1} = zetaCentr;
        theta_cell{exportCounter,1} = theta(h_ad.val - zetaCentr);
        flux_cell{exportCounter,1} = q(h_ad.val,h_m0);
        t(exportCounter,1) = timeCum;
        iterations(exportCounter,1) = iter-1;
        cpuTime(exportCounter,1) = tf;
        exportCounter = exportCounter + 1;   
    end

end

%% Storing results in the sol structure
sol.time = t;           sol.h = h_cell;             sol.psi = psi_cell;
sol.zeta = zeta_cell;   sol.theta = theta_cell;     sol.flux = flux_cell;
sol.iter = iterations;  sol.residual = residual;    sol.cpuTime = cpuTime;

fprintf('\n sol: \n'); disp(sol);  % printing solutions in the command window

%% Water content as a function of time
figure(3);
for ii=1:printLevels
    plotCellData(G,sol.theta{ii},sol.theta{ii} > theta(-990));
    xlabel('x-axis [cm]'), ylabel('y-axis [cm]'), zlabel('Depth [cm]');
    xlim([0 Lx]); ylim([0 Ly]); zlim([0 Lz]);
    view([-56 -51]); shading faceted; camproj perspective;
    set(gca, 'ZDir', 'reverse'); box on; pbaspect([1 1 2.5]);
    t_h1=handle(text('Units','normalized',...    % Use [%] of the axis length
                     'Position',[1,0.9,1],...    % Position of text   
                     'EdgeColor','w'));          % Textbox 
    t_h1.String=sprintf('Time: %2.2f [h]',sol.time(ii)/3600); % print "time" counter
    title('Water Content > 0.11'); pause(0.5);
    colorbar; caxis([theta(psiB) theta(psiT)]);
    t_h1.String =[];
end