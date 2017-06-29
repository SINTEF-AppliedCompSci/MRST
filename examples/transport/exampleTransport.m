clear;

mrstModule add ad-core ad-props ad-blackoil geochemistry mrst-gui


%% Define the grid

G = cartGrid([100, 1, 1], [10, 1, 1]);
G = computeGeometry(G);

% plotGrid(G), view(3), axis tight

%% Define the rock

rock.perm = 1*darcy*ones(G.cells.num, 1);
rock.poro = ones(G.cells.num, 1);

%% Define the fluid
pRef = 1*barsa;
fluid = initSimpleADIFluid('phases', 'W', 'mu', 1*centi*poise, 'rho', ...
                           1000*kilogram/meter^3, 'c', 1e-10, 'cR', 4e-10, ...
                           'pRef', pRef);
fluid = initSimpleADIFluid('phases', 'W', 'mu', 1*centi*poise, 'rho', ...
                           1000*kilogram/meter^3, 'c', 0, 'cR', 0, 'pRef', ...
                           pRef);

%% Define the chemistry
MasterCompNames = {'O', 'H', 'Na*', 'Cl*'};
CompNames = {'H+*', 'OH-', 'Na+', 'Cl-', 'NaCl', 'H2O*'};
Reactions ={'H2O  <-> H+  + OH- ',      10^-14*mol/litre, ...
            'NaCl <-> Na+ + Cl-',       10^1*mol/litre};
                 
% instantiate chemical model
chemModel = ChemicalModel(CompNames, MasterCompNames, Reactions);
chemModel.printChemicalSystem;

% initial chemistry
inputConstraints = [1e-3 1e-3 1e-9 1]*mol/litre;
[initchemstate, initreport]= chemModel.initState(inputConstraints);

% injected chemistry
inputConstraints = [1e-2 1e-2 1e-3 1]*mol/litre;
[injchemstate, injreport] = chemModel.initState(inputConstraints);


%% Define the initial state

nc = G.cells.num;
initState.components          = repmat(initchemstate.components, nc, 1);
initState.masterComponents    = repmat(initchemstate.masterComponents, nc, 1);
initState.logcomponents       = repmat(initchemstate.logcomponents, nc, 1);
initState.logmasterComponents = repmat(initchemstate.logmasterComponents, nc, 1);

initState.pressure         = linspace(pRef, 0, nc + 1)';
initState.pressure(end)    = [];


%% Define the boundary conditions

src                  = [];
src                  = addSource(src, [1], [1/100].*meter^3/day, 'sat', 1);
src.masterComponents = injchemstate.masterComponents; 
src.logmasterComponents = injchemstate.logmasterComponents; 

bc                     = [];
bc                     = pside(bc, G, 'east', 0*barsa, 'sat', 1);
bc.masterComponents    = [initchemstate.masterComponents];        % (will not used if outflow)
bc.logmasterComponents = [initchemstate.logmasterComponents];     % (will not used if outflow)

%% Define the model

model = ChemicalTransportSplitLogModel(G, rock, fluid, chemModel);

%% Define the schedule

schedule.step.val     =  [1*day*ones(5, 1); 10*day*ones(40, 1)];
schedule.step.control = ones(numel(schedule.step.val), 1);
schedule.control = struct('bc', bc, 'src', src, 'W', []);

% We need to figure out how the control is turned into a drivingForces,

%% Run the simulation

% timestepselector = IterationCountTimeStepSelector('targetIterationCount', 10);
% solver = NonLinearSolver('timeStepSelector', timestepselector);

solver = NonLinearSolver('minIterations', 1);

[~, states, schedulereport] = simulateScheduleAD(initState, model, schedule, ...
                                                 'NonLinearSolver', solver);


statesNewvar = change_units(states, litre/mol);

figure
plotToolbar(G, statesNewvar, 'plot1d', true, 'log10', true);

