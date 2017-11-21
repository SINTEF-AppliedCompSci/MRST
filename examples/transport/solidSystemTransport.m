close all;
clear;

mrstModule add ad-core ad-props ad-blackoil geochemistry

mrstVerbose off

%% Define the grid

G = cartGrid([100, 1, 1], [10, 1, 1]);
G = computeGeometry(G);
nc = G.cells.num;

%% Define the rock

rock.perm = 1*darcy*ones(G.cells.num, 1);
rock.poro = ones(G.cells.num, 1);

%% Define the fluid
pRef = 0*barsa;
fluid = initSimpleADIFluid('phases', 'W', 'mu', 1*centi*poise, 'rho', ...
                           1000*kilogram/meter^3, 'c', 1e-10, 'cR', 4e-10, ...
                           'pRef', pRef);

%% Define the chemistry

elements = {'O', 'H', 'Na*', 'Cl*', 'Ca*', 'C*'};

species = {'H+*', 'OH-', 'Na+', 'Cl-', 'NaOH', 'H2O*',...
             'Ca+2', 'CO3-2', 'HCO3-', 'CO2',...
             'CaCO3(s)'};

reactions ={'H2O  = H+  + OH- ',            10^-14*mol/litre, ...
            'NaOH = Na+ + OH-',             10^10*mol/litre,...
            'CaCO3(s) = CO3-2 + Ca+2',      10^-8.48*(mol/litre)^2,...
            'CO3-2 + H+ = HCO3-',           10^10.329/(mol/litre),...
            'CO3-2 + 2*H+ = CO2 + H2O',     10^16.681/(mol/litre)};

% instantiate chemical model object
chemModel = ChemicalModel(elements, species, reactions);

% print the chemical model to the screen
chemModel.printChemicalSystem;

% initial chemistry
Nai = 1e-1;
Cli = Nai;
Cai = 1e-2;
Ci = 1e-3;
Hi = 1e-9;
H2Oi = 1;

inputConstraints = [Nai Cli Cai Ci Hi H2Oi]*mol/litre;
[initchemstate, initreport]= chemModel.initState(repmat(inputConstraints, nc,1), 'charge', 'Cl');

% injected chemistry
Naf = 1e-1;
Clf = Naf;
Caf = 1e-3;
Cf = 1e-3;
Hf = 1e-9;
H2Of = 1;

inputConstraints = [Naf Clf Caf Cf Hf H2Of]*mol/litre;
[injchemstate, injreport] = chemModel.initState(inputConstraints, 'charge', 'Cl');

%% Define the initial state

initState = initchemstate;
initState.pressure = pRef*ones(nc,1);

%% Define the transport model
model = ChemicalTransportModel(G, rock, fluid, chemModel);

% plot the species and element distribution at the end of each time step
model.plotFinal = false;
model.plotIter = false;

% plotting during the simulation slows the solver down quite a bit

%% Define the boundary conditions

% use model.fluidMat to pull the fluid concentrations from the injected
% state
injfluidpart = injchemstate.species*model.fluidMat';
initfluidpart = initchemstate.species(end,:)*model.fluidMat';

pv = poreVolume(G,rock);

% define source term at cell 1 for inflow boundary condition
src                	= [];
src               	= addSource(src, 1, pv(1)/day, 'sat', 1);

% give the fluid concentration at the inlet
src.elements        = injfluidpart;
src.logElements   	= log(injfluidpart);

% give dirchlet boundary condition at outlet
bc                  = [];
bc                  = pside(bc, G, 'east', 0*barsa, 'sat', 1);
bc.elements         = initfluidpart;        % (will not be used if outflow)
bc.logElements      = log(initfluidpart);  % (will not be used if outflow)


%% Define the schedule

% ten time steps of 0.01 days followed by 100 steps of 1 day
schedule.step.val = [0.01*day*ones(10, 1); 1*day*ones(100, 1);];
schedule.step.control = ones(numel(schedule.step.val), 1);
schedule.control = struct('bc', bc, 'src', src, 'W', []);


%% Run the simulation

[~, states, scheduleReport] = simulateScheduleAD(initState, model, schedule);

plotToolbar(G, states)