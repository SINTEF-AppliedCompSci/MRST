mrstModule add ad-core ad-eor ad-blackoil ad-props blackoil-sequential matlab_bgl

gravity reset on

n = 10;
G = computeGeometry(cartGrid([n,1,1]));
rock = makeRock(G, 1, 1);

fluid = initSimpleADIFluid('n'     , [2, 2, 2], ...
                           'rho'   , [1000, 800, 100]*kilogram/meter^3, ...
                           'phases', 'WOG', ...
                           'mu'    , [1, 10, 0.1]*centi*poise);

sres = 0;      
fluid.krO = coreyPhaseRelpermAD(2, sres, 1, sres);
fluid.sOres = sres;
fluid.sGres = 0;
fluid.mixPar = 0;

model = OilWaterSolventModel2(G, rock, fluid);

solventPeriod = 0.2;

T = 1*year;

Tsol = solventPeriod*T;
pv = poreVolume(G, rock);
injRate = 1*sum(pv)/Tsol;

%%

W = [];

W = addWell(W, G, rock, 1, 'comp_i', [0, 0, 1], 'type', 'rate', 'val', injRate);
W = addWell(W, G, rock, G.cells.num, 'comp_i', [0, 0, 1], 'type', 'bhp' , 'val', 0      );

state0 = initResSol(G, 100*barsa, [0 1 0]);
state0.wellSol = initWellSolAD(W, model, state0);

nStep = 100;
% dT = repmat(Tsol/nStep, nStep, 1);
dT = rampupTimesteps(Tsol, Tsol/nStep);
schedule = simpleSchedule(dT, 'W', W);

nls = NonLinearSolver('useLinesearch', false);
[ws, statesS, reports] = simulateScheduleAD(state0, model, schedule, 'nonlinearsolver', nls);


%%

Twat = (1-solventPeriod)*T;
injRate = 0.2*sum(pv)/Twat;

W = [];
W = addWell(W, G, rock, 1, 'comp_i', [1, 0, 0], 'type', 'rate', 'val', injRate);
W = addWell(W, G, rock, G.cells.num, 'comp_i', [0, 0, 1], 'type', 'bhp' , 'val', 0      );

nStep = 100;
dT = repmat(Twat/nStep, nStep, 1);
schedule = simpleSchedule(dT, 'W', W);

state0 = statesS{end};

[ws, statesW, reports] = simulateScheduleAD(state0, model, schedule);

%%

figure(1); clf
states = cat(1,statesS, statesW);
plotToolbar(G, states, 'plot1d', true)

%%

W = [];

injRate = 0.2*sum(pv)/Twat;

W = addWell(W, G, rock, 1, 'comp_i', [1, 0, 0], 'type', 'rate', 'val', injRate);
W = addWell(W, G, rock, G.cells.num, 'comp_i', [0, 0, 1], 'type', 'bhp' , 'val', 0      );

state0 = initResSol(G, 100*barsa, [0 1 0]);
state0.wellSol = initWellSolAD(W, model, state0);

nStep = 100;
dT = repmat(Twat/nStep, nStep, 1);
schedule = simpleSchedule(dT, 'W', W);

[wsWW, statesWW, reports] = simulateScheduleAD(state0, model, schedule);

%%

figure(2); clf;
plotToolbar(G, statesWW, 'plot1d', true)
