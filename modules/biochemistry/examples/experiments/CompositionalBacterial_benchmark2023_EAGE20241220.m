%% MRST Simulation for Hydrogen Storage with Bacterial Growth Model
% Description: This script uses MRST to model gas injection into a 3D porous medium,
% incorporating compositional fluid properties, bacterial mono modal.
% We consider a liquid phase (W) and a gas (G) phase, 4 components 
% ('H2O','H2','CO2','CH4') and The microbial activity of methanogenic archaea.
%This test case comes from a Benchmark in EAGE 2023

% Clear workspace and initialize MRST modules
clear; clc;
mrstModule add biochemistry compositional ad-blackoil ad-core ad-props mrst-gui
gravity reset on

%% ============Grid and Rock Properties=====================
% Define grid dimensions and physical dimensions
%[nx, ny, nz] = deal(61,61,10);  % Grid cells in x, y, z directions
[nx, ny, nz] = deal(31,31,8);  % Grid cells in x, y, z directions
[Lx,Ly,Lz] = deal(1525,1525,50);         % Physical dimensions in meters
dims = [nx, ny, nz];
pdims = [Lx, Ly, Lz];


% Create grid and shift vertically by reservoir depth
G = cartGrid(dims, pdims);
depth_res = 3368;                 % Reservoir depth in meters
G.nodes.coords(:, 3) = G.nodes.coords(:, 3) + depth_res;
G = computeGeometry(G);

% Define rock properties
K=[100, 100, 10].*milli*darcy;
rock = makeRock(G, K, 0.2);  % Default permeability and porosity


%% Fluid Properties Initialization
% Define compositional fluid model (with CoolProp library support)
compFluid = TableCompositionalMixture({'Water', 'Hydrogen', 'CarbonDioxide', 'Methane'}, ...
                                      {'H2O', 'H2', 'CO2', 'C1'});

% Fluid density and viscosity (kg/m^3 and cP)
[rhow, rhog] = deal(996.52 * kilogram / meter^3, 38.3823 * kilogram / meter^3);
[viscow, viscog] = deal(0.65403 * centi * poise, 0.01123 * centi * poise);

% Compressibility (per bar)
[cfw, cfg] = deal(0, 9.031e-3 / barsa);

% Relative permeability and initial saturations
[srw, src] = deal(0.0, 0.0);
P0=100 * barsa;
fluid = initSimpleADIFluid('phases', 'WG', 'mu', [viscow, viscog], ...
                           'rho', [rhow, rhog], 'pRef', P0, ...
                           'c', [cfw, cfg], 'n', [2, 2], 'smin', [srw, src]);

% Capillary pressure function
Pe = 0.1 * barsa;
pcWG = @(sw) Pe * sw.^(-1/2);
fluid.pcWG = @(sg) pcWG(max((1 - sg - srw) / (1 - srw), 1e-5));

%% Simulation Parameters
% Set total time, pore volume, and injection rate
niter=60;
TotalTime = niter*day;% 105*day; 
rate = 1e6*meter^3/day; 
%pv = sum(poreVolume(G, rock)) / TotalTime;
%rate = 10*pv; %100 * pv;



%% Time Stepping and Schedule
% Define schedule and solver
nls = NonLinearSolver('useRelaxation', true);
deltaT =rampupTimesteps(TotalTime, 1*day, 0);
schedule = simpleSchedule(deltaT);
nj1=30;
schedule.step.control(1:nj1)=1;
schedule.step.control(nj1+1:end)=2;

%% Wells and Boundary Conditions
% Initialize wells
W1 = [];
W2 = [];
tmp = cell(2,1);
n1=floor(0.5*nx)+1; n2=floor(0.5*nx)+1;
schedule.control = struct('W',tmp);

% Injection well parameters
W1 = verticalWell(W1, G, rock, n1, n2, nz, 'comp_i', [0, 1], 'Radius', 0.5, ...
                 'name', 'Injector', 'type', 'rate', 'Val', rate, 'sign', 1);
W1.components = [0.0, 0.95,  0.05, 0.0];  % H2-rich injection   {'H2O', 'H2', 'CO2', 'C1'});
W2 = verticalWell(W2, G, rock, n1, n2, nz, 'comp_i', [0, 1], 'Radius', 0.5, ...
                 'name', 'Rest', 'type', 'rate', 'Val', 0.0, 'sign', 1);
W2.components = [0.0, 0.95,  0.05, 0.0];  % H2-rich injection   {'H2O', 'H2', 'CO2', 'C1'});


schedule.control(1).W = W1;
schedule.control(2).W = W2;


%% Model Setup: Compositional Model with Bacterial Growth
arg = {G, rock, fluid, compFluid, 'water', true, 'oil', false, 'gas', true, ...
       'bacteriamodel', false, 'bDiffusionEffect', false, ...
       'moleculardiffusion',false,'liquidPhase', 'W', 'vaporPhase', 'G'};
model = BiochemistryModel(arg{:});
model.outputFluxes = false;
eosname='pr';% 'sw';
model.EOSModel = SoreideWhitsonEquationOfStateModel(G, compFluid,eosname);
model.EOSModel.msalt=0;

%% Initial Conditions
% Temperature and initial saturations
T0 = 313.15;                % Initial temperature (K)
s0 = [0.2, 0.8];           % Initial saturations (Sw,Sg)
z0 = [0.2, 0.0, 0.0, 0.8];  % Initial composition: H2O, H2, CO2, CH4
Phydro0=rhow*norm(gravity).*G.cells.centroids(:,3);
% Initialize state with bacterial concentration
nbact0 = 10^6;
state0 = initCompositionalStateBacteria(model, Phydro0, T0, s0, z0, nbact0);


%% Run simulation
[~, states, report] = simulateScheduleAD(state0, model, schedule,'nonlinearsolver', nls);


%% Plotting Results
namecp = model.EOSModel.getComponentNames();
indH2=find(strcmp(namecp,'H2'));
indCO2= find(strcmp(namecp,'CO2'));
nT=numel(states);
zH2=zeros(nT,1);
zCO2= zeros(nT,1);
for i = 1:nT
    zH2(i) = max(states{i}.components(:,indH2));
    zCO2(i) = max(states{i}.components(:,indCO2));
end

for i = 1:nT
    figure(1); clf; hold on
    plot(1:nT,zH2,'b',1:nT,zCO2,'k-')
end
for i = 1:nT
    figure(2); clf; hold on
    plot(1:nT,zCO2,'k-')
end
