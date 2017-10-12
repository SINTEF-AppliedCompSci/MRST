

%% Load part of the SPE10 model

mrstModule add spe10

% We choose 5 layers from the Tarbert formation
[G, W, rock] = SPE10_setup((5:9)');


%% Load fluid data with polymer properties

fn    = 'POLYMER.DATA';
deck  = readEclipseDeck(fn);
deck  = convertDeckUnits(deck);
fluid = initDeckADIFluid(deck);
clear deck


%% Create models

modelFI = OilWaterPolymerModel(G, rock, fluid);
modelSQ = getSequentialModelFromFI(modelFI);




%% Wells

W = addWell([], G, rock, find(ijk{1}==10&ijk{2}==25), 'type', 'rate', ...
   'val', 200*meter^3/day, 'radius', 0.1, 'sign', 1, 'name', 'INJE', ...
   'comp_i', [1 0]);
W = addWell(W, G, rock, find(ijk{1}==90&ijk{2}==25), 'type', 'bhp', ...
   'val', 250*barsa, 'radius', 0.1, 'sign', -1, 'name', 'PROD', ...
   'comp_i', [0 1]);


%% Initial state

state0 = initState(G, W, 250*barsa, [swir 1-swir]);

%% Schedule


%state0.wellSol = initWellSolAD(W, model, state0);


% Schedule
schedule = simpleSchedule(1*day.*ones(10,1), 'wells', W);


