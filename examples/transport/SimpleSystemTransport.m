close all;
clear;

mrstModule add ad-core ad-props ad-blackoil geochemistry mrst-gui

%% Define the grid
sideLength = 50;
G = cartGrid([sideLength, sideLength, 1], [10, 10, 1]);
G = computeGeometry(G);
nc = G.cells.num;

plotGrid(G), view(3), axis tight
outInd = [G.faces.num-G.cells.num-sideLength+1; G.faces.num-G.cells.num-sideLength^2+sideLength];

plotFace(G, outInd, 'r');

%% Define the rock

rock.perm = 1*darcy*ones(G.cells.num, 1);
rock.poro = 0.5*ones(G.cells.num, 1);

%% Define the fluid
pRef = 0*barsa;
fluid = initSimpleADIFluid('phases', 'W', 'mu', 1*centi*poise, 'rho', ...
                           1000*kilogram/meter^3, 'c', 1e-10, 'cR', 4e-10, ...
                           'pRef', pRef);

%% Define the chemistry

elements = {'O', 'H', 'Na*','Cl*'};

species = {'H+*', 'OH-', 'Na+', 'H2O*', 'NaCl','Cl-'};

reactions ={'H2O  = H+  + OH- ',      10^-14*mol/litre, ... 
            'NaCl = Na+ + Cl-',       10^1*mol/litre};


% instantiate chemical model
chemModel = ChemicalModel(elements, species, reactions);

% initial chemistry
Nai = 1e-3;
Cli = Nai;
Hi = 1e-7;
H2Oi = 1;
initial = [Nai Cli Hi H2Oi]*mol/litre;
[initChemState, initreport]= chemModel.initState(repmat(initial, nc,1), 'charge', 'Cl');

% injected chemistry
Naf = 1e-3;
Clf = Naf;
Hf = 1e-10;
H2Of = 1;
injected = [Naf Clf Hf H2Of]*mol/litre;
[injChemState, injreport]= chemModel.initState(injected, 'charge', 'Cl');

%% Define the initial state

initChemState.pressure          = pRef*ones(nc,1);

%% Define the model
model = ChemicalTransportModel(G, rock, fluid, chemModel);
model.plotIter = false;

%% Define the boundary conditions

pv = poreVolume(G,rock);


src                	= [];
src               	= addSource(src, [1; nc], pv(1:2)/day, 'sat', 1);
src.elements        = injChemState.elements(end,:);
src.logElements     = injChemState.logElements(end,:);

bc                  = [];
outInd = [G.faces.num-G.cells.num-sideLength+1; G.faces.num-G.cells.num-sideLength^2+sideLength];
bc                  = addBC(bc, outInd, 'pressure', [0; 0]*barsa, 'sat', 1);
bc.elements         = initChemState.elements(end,:);        % (will not used if outflow)
bc.logElements      = initChemState.logElements(end,:);  % (will not used if outflow)


%% Define the schedule

schedule.step.val = [0.01*day*ones(5, 1); 0.1*day*ones(5,1); 1*day*ones(5, 1); 5*day*ones(10, 1)];
schedule.step.control = ones(numel(schedule.step.val), 1);
schedule.control = struct('bc', bc, 'src', src, 'W', []);


%% Run the simulation

[~, states, scheduleReport] = simulateScheduleAD(initChemState, model, schedule);

% visualize
plotToolbar(G, states)