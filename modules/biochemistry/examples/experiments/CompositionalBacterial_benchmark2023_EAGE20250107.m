%% MRST Simulation for Hydrogen Storage with Bacterial Growth Model
% Description: This script uses MRST to model gas injection into a 3D porous medium,
% incorporating compositional fluid properties, bacterial mono modal.
% We consider a liquid phase (W) and a gas (G) phase, 4 components 
% ('H2O','H2','CO2','CH4') and The microbial activity of methanogenic archaea.
%This test case comes from a Benchmark in EAGE 2023
%TO DO: ERROR BUG WITH BACTERIAL=FALSE AND PRODUCTION WELL. DIVERGES WITH
%BACTERIA
% Clear workspace and initialize MRST modules
clear; clc;
%mrstModule add biochemistry compositional ad-blackoil ad-core ad-props mrst-gui
mrstModule add biochemistry compositional ad-blackoil ad-core ad-props mrst-gui
gravity reset on
biochemistrymodel=true;
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
niter=105;% 22;
TotalTime = niter*day;
rate = 1e6*meter^3/day; 



%% Time Stepping and Schedule
% Define schedule and solver
nls = NonLinearSolver('useRelaxation', true);
deltaT =rampupTimesteps(TotalTime, 1*day, 0);
schedule = simpleSchedule(deltaT);
nj1=30;nj2=60;nj3=90;
%nj1=5;nj2=10;nj3=15;
schedule.step.control(1:nj1)=1;
schedule.step.control(nj1+1:nj2)=2;
schedule.step.control(nj2+1:nj3)=3;
schedule.step.control(nj3+1:end)=4;

%% Wells and Boundary Conditions
% Initialize wells
W1 = [];
W2 = [];
W3 = [];
W4 = [];
tmp = cell(4,1);
n1=floor(0.5*nx)+1; n2=floor(0.5*nx)+1;
schedule.control = struct('W',tmp);

% Injection well parameters
W1 = verticalWell(W1, G, rock, n1, n2, nz, 'comp_i', [0, 1], 'Radius', 0.5, ...
                 'name', 'Injector', 'type', 'rate', 'Val', rate, 'sign', 1);
W1(1).components = [0.0, 0.95,  0.05, 0.0];  % H2-rich injection   {'H2O', 'H2', 'CO2', 'C1'});

%Idle period
W2 = verticalWell(W2, G, rock, n1, n2, nz, 'comp_i', [0, 1], 'Radius', 0.5, ...
                 'name', 'Rest', 'type', 'rate', 'Val', 0.0, 'sign', 1);
W2(1).components = [0.0, 0.95,  0.05, 0.0];  % rest period

%production
W3 = verticalWell(W3, G, rock, n1, n2, nz, 'comp_i', [0, 1], 'Radius', 0.5, ...
                 'name', 'Prod', 'type', 'rate', 'Val', -rate, 'sign', -1);
W3(1).components = [0.0, 0.95,  0.05, 0.0];  %production
W3.lims.bhp= P0;

%Idle period
W4 = verticalWell(W4, G, rock, n1, n2, nz, 'comp_i', [0, 1], 'Radius', 0.5, ...
                 'name', 'Rest', 'type', 'rate', 'Val', 0.0, 'sign', 1);
W4(1).components = [0.0, 0.95,  0.05, 0.0];  % rest period


schedule.control(1).W = W1;
schedule.control(2).W = W2;
schedule.control(3).W = W3;
schedule.control(4).W = W4;


%% Model Setup: Compositional Model with Bacterial Growth
if biochemistrymodel
    arg = {G, rock, fluid, compFluid, 'water', true, 'oil', false, 'gas', true, ...
       'bacteriamodel', true, 'bDiffusionEffect', false, ...
       'moleculardiffusion',false,'liquidPhase', 'W', 'vaporPhase', 'G'};
    model = BiochemistryModel(arg{:});
    model.outputFluxes = false;
    eosname='pr';% 'sw';
    model.EOSModel = SoreideWhitsonEquationOfStateModel(G, compFluid,eosname);
    model.EOSModel.msalt=0;
else
    arg = {G, rock, fluid, compFluid, 'water', true, 'oil', false, 'gas', true, ...
      'liquidPhase', 'W', 'vaporPhase', 'G'};
    model = GenericOverallCompositionModel(arg{:});
end



%% Initial Conditions
% Temperature and initial saturations
T0 = 313.15;                % Initial temperature (K)
s0 = [0.2, 0.8];           % Initial saturations (Sw,Sg)
z0 = [0.2, 0.0, 0.0, 0.8];  % Initial composition: H2O, H2, CO2, CH4
Phydro0=rhow*norm(gravity).*G.cells.centroids(:,3);
% Initialize state with bacterial concentration

if biochemistrymodel
    if model.bacteriamodel
        nbact0 = 10^6;
        state0 = initCompositionalStateBacteria(model, P0, T0, s0, z0, nbact0);
    else
        state0 = initCompositionalState(model, P0, T0, s0, z0);
    end

else
    state0 = initCompositionalState(model, P0, T0, s0, z0);
end


%% Run simulation
%mrstModule add mpfa
%model_mpfa = setMPFADiscretization(model);
[wellSols, states, report] = simulateScheduleAD(state0, model, schedule,'nonlinearsolver', nls);


%% Plotting Results
%plotWellSols(wellSols, cumsum(schedule.step.val))
namecp = model.EOSModel.getComponentNames();
indH2=find(strcmp(namecp,'H2'));
indCO2= find(strcmp(namecp,'CO2'));
indCH4= find(strcmp(namecp,'C1'));
nT=numel(states);
xH2=zeros(nT,1);
yH2=zeros(nT,1);
yCO2= zeros(nT,1);
yCH4= zeros(nT,1);
for i = 1:nT
    xH2(i)=max(states{i}.x(:,indH2));
    yH2(i)=max(states{i}.y(:,indH2));
    yCO2(i)=max(states{i}.y(:,indCO2));
    yCH4(i)=max(states{i}.y(:,indCH4));
end

for i = 1:nT
    figure(1); clf; 
    plot(1:nT,yH2,'b')
    hold on
    plot(1:nT,yCO2,'k-')
    hold on
    plot(1:nT,yCH4,'r')
end
title('Molar fractions')
xlabel('Time (days)')
ylabel('molar fraction')
legend('yH2','yCO2','yCH4')
