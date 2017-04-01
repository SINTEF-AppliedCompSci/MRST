mrstModule add ad-core ad-eor ad-blackoil ad-props blackoil-sequential matlab_bgl

G = computeGeometry(cartGrid([10,1,1]));
rock = makeRock(G, 1, 1);

fluid = initSimpleADIFluid('n', [2, 2, 2], ...
                           'rho', [1, 1, 1], ...
                           'phases', 'WOG', ...
                           'mu', [1, 1, 1]*centi*poise);

fluid = addSolventProperties(fluid);

                       
T = 1*year;
pv = poreVolume(G, rock);
injRate = sum(pv)/T;
W = [];
W = addWell(W, G, rock, 1, 'comp_i', [0, 0, 1], 'type', 'rate', 'val', injRate);
W = addWell(W, G, rock, G.cells.num, 'comp_i', [0, 0, 1], 'type', 'bhp', 'val', 0);

state0 = initState(G, W, 0, [1,0,0]);

model = OilWaterSolventModel(G, rock, fluid);

nStep = 10;
dT = repmat(T/nStep, nStep, 1);

schedule = simpleSchedule(dT, 'W', W);

[ws, states, reports] = simulateScheduleAD(state0, model, schedule);