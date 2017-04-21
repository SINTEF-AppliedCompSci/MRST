mrstModule add ad-core ad-eor ad-blackoil ad-props blackoil-sequential matlab_bgl

%%------------------------------------------------------------------------------------

gravity reset off

n = 20;
G = computeGeometry(cartGrid([n,n,1]));
rock = makeRock(G, 1, 1);

fluid = initSimpleADIFluid('n'     , [2, 2, 2], ...
                           'rho'   , [1000, 800, 0.1]*kilogram/meter^3, ...
                           'phases', 'WOG', ...
                           'mu'    , [100, 1000, 0.01]*centi*poise);

sres = 0;      
fluid.krO = coreyPhaseRelpermAD(2, sres, 1, sres);
fluid.sOres = sres;
fluid.sGres = 0;
fluid.mixPar = 1;

model = OilWaterSolventModel2(G, rock, fluid);

T = 1*year;

pv = poreVolume(G, rock);
injRate = 1*sum(pv)/T;

%%------------------------------------------------------------------------------------

W = [];

W = addWell(W, G, rock, 1, 'comp_i', [0, 0, 1], 'type', 'rate', 'val', injRate);
W = addWell(W, G, rock, G.cells.num, 'comp_i', [0, 0, 1], 'type', 'bhp' , 'val', 0      );

state0 = initResSol(G, 100*barsa, [0 1 0]);
state0.wellSol = initWellSolAD(W, model, state0);

nStep = 100;
dT = repmat(T/nStep, nStep, 1);
schedule = simpleSchedule(dT, 'W', W);

[ws, states, reports] = simulateScheduleAD(state0, model, schedule);


%%------------------------------------------------------------------------------------

figure(1); clf
plotToolbar(G, states)
