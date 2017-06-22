mrstModule add ad-core ad-eor ad-blackoil ad-props blackoil-sequential matlab_bgl

gravity reset on

n = 200;
G = computeGeometry(cartGrid([n,1,1]));
rock = makeRock(G, 100*milli*darcy, 1);

fluid = initSimpleADIFluid('n'     , [2, 2, 2], ...
                           'rho'   , [1000, 800, 100]*kilogram/meter^3, ...
                           'phases', 'WOG', ...
                           'mu'    , [1, 10, 1]*centi*poise);
                       
%'c'     , [1e-7, 1e-6, 1e-6]/barsa, ...
    
sres = 0;      
fluid.sOres = sres;
fluid.sGres = 0;
fluid.mixPar = 1;

model = MiscibleOilWaterSolventModel(G, rock, fluid);

T = 1*year;
pv = poreVolume(G, rock);
injRate = 2*sum(pv)/T;

nStep = 100;
dT = rampupTimesteps(T, T/nStep);

%%

W = [];
W = addWell(W, G, rock, 1, 'comp_i', [0.5, 0, 0.5], 'type', 'rate', 'val', injRate);
W = addWell(W, G, rock, G.cells.num, 'comp_i', [0, 0, 1], 'type', 'bhp' , 'val', 0      );
schedule = simpleSchedule(dT, 'W', W);

state0 = initResSol(G, 100*barsa, [0 1 0]);
state0.wellSol = initWellSolAD(W, model, state0);

[wsS, statesS, reportsS] = simulateScheduleAD(state0, model, schedule);

%%

mrstModule add mrst-gui

figure(1); clf
plotToolbar(G, statesS, 'plot1d', true)
ylim([0 1]);

%%

W = [];
W = addWell(W, G, rock, 1, 'comp_i', [1, 0, 0], 'type', 'rate', 'val', injRate);
W = addWell(W, G, rock, G.cells.num, 'comp_i', [0, 0, 1], 'type', 'bhp' , 'val', 0      );
schedule = simpleSchedule(dT, 'W', W);

state0 = initResSol(G, 100*barsa, [0 1 0]);
state0.wellSol = initWellSolAD(W, model, state0);

[wsW, statesW, reportsW] = simulateScheduleAD(state0, model, schedule);

%%

figure(2); clf
plotToolbar(G, statesW, 'plot1d', true)
ylim([0 1]);
