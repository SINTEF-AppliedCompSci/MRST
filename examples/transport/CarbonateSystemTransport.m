close all;
clear;

mrstModule add ad-core ad-props ad-blackoil geochemistry mrst-gui


%% Define the grid

G = cartGrid([10, 1, 1]);
G = computeGeometry(G);

plotGrid(G), view(3), axis tight

%% Define the rock

rock.perm = 1*darcy*ones(G.cells.num, 1);
rock.poro = ones(G.cells.num, 1);

%% Define the fluid
pRef = 0*barsa;
fluid = initSimpleADIFluid('phases', 'W', 'mu', 1*centi*poise, 'rho', ...
                           1000*kilogram/meter^3, 'c', 1e-10, 'cR', 4e-10, ...
                           'pRef', pRef);

%% Define the chemistry

chemModel = callCarbonateChemicalSystem();

% initial chemistry
inputConstraints = [1 1e-2 1e-7 1e-3]*mol/litre;
[initchemstate, initreport]= chemModel.initState(inputConstraints);

% injected chemistry
inputConstraints = [1 1e-2 1e-9 1e-3]*mol/litre;
[injchemstate, injreport] = chemModel.initState(inputConstraints);

%% Define the initial state

nc = G.cells.num;
initState.components       = repmat(initchemstate.components, nc, 1);
initState.masterComponents = repmat(initchemstate.masterComponents, nc, 1);

% initState.components([1],:)         = repmat(injchemstate.components,1,1);
% initState.masterComponents([1],:)	= repmat(injchemstate.masterComponents,1,1);

initState.pressure         = pRef*ones(nc,1);


%% Define the boundary conditions


src                  = [];
src                  = addSource(src, [1], [1].*meter^3/day, 'sat', 1);
src.masterComponents = injchemstate.masterComponents; % concentration composition


bc                  = [];
bc                  = pside(bc, G, 'east', 0*barsa, 'sat', 1);
bc.masterComponents= [initchemstate.masterComponents];                                                 % (will not used if outflow)

W = []; % this is a hack: the field seems to be needed for the moment.

%% Define the model

model = ChemicalTransportModel(G, rock, fluid, chemModel);

%% Define the schedule

schedule.step.val     =  1*day*ones(30, 1);
schedule.step.control = ones(numel(schedule.step.val), 1);
schedule.control = struct('bc', bc, 'src', src, 'W', W);

% We need to figure out how the control is turned into a drivingForces,

%% Run the simulation

[~, states] = simulateScheduleAD(initState, model, schedule);

[ states ] = change_units( states, litre/mol );

plotToolbar(G, states,'plot1d', true, 'log10', true);