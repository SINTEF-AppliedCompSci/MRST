mrstModule add mrst-solvent spe10

%% Use setupSPE10_AD to Fetch the SPE10 model
% We pick up only one layer 
%
layers = 35;
[~, model, ~] = setupSPE10_AD('layers', layers);
% We recover the grid and rock properties from the model
G = model.G;
rock = model.rock;

%%

gravity reset on

fluid = initSimpleADIFluid('n'     , [2, 2, 2], ...
                           'rho'   , [1000, 800, 100]*kilogram/meter^3, ...
                           'phases', 'WOG', ...
                           'mu'    , [1, 10, 2]*centi*poise);
%'c'     , [1e-7, 1e-6, 1e-6]/barsa, ...
                       
fluid = addSolventProperties(fluid, 'n', 2, ...
                                    'rho', 100*kilogram/meter^3, ...
                                    'mixPar', 2/3, ...
                                    'mu'    , 1*centi*poise);
%'c', 1e-6/barsa, ...                                

model = FourPhaseSolventModel(G, rock, fluid);

T = 1*year;
pv = poreVolume(G, rock);
injRate = 1*sum(pv)/T;

nStep = 200;
dT = T/nStep;
%%

WA = verticalWell([], G, rock, 1, 1, [], ...
                 'comp_i', [0, 0, 0, 1], ...
                 'type'  , 'rate', ...
                 'val'   , injRate);
     
WA = verticalWell(WA, G, rock, 60, 220, [], ...
                 'comp_i', [0, 0, 0, 1], ...
                 'type', 'bhp', ...
                 'val', 0      );
schedule = simpleSchedule(dT, 'W', WA);

WB = verticalWell([], G, rock, 1, 1, [], ...
                 'comp_i', [1, 0, 0, 0], ...
                 'type'  , 'rate', ...
                 'val'   , injRate);

WB = verticalWell(WB, G, rock, 60, 220, [], ...
                 'comp_i', [0, 0, 0, 1], ...
                 'type'  , 'bhp', ...
                 'val'   , 0);
             
control(1).W = WA;
control(2).W = WB;

injStart = 0;
injStop = 0.25;

dT_A = rampupTimesteps(injStop*T, dT);
dT_B = rampupTimesteps((1 - injStop)*T, dT);

step.val = [dT_A; dT_B];
step.control = [1*ones(numel(dT_A),1); 2*ones(numel(dT_B),1)];
schedule.control = control;
schedule.step = step;

state0 = initResSol(G, 100*barsa, [0 1 0 0]);
state0.wellSol = initWellSolAD(WA, model, state0);

%%

[ws, states, reports] = simulateScheduleAD(state0, model, schedule);


%%

plotToolbar(G, states)