mrstModule add ad-core ad-eor ad-blackoil ad-props blackoil-sequential matlab_bgl

gravity reset on

n = 10;
G = computeGeometry(cartGrid([n,1,1]));
rock = makeRock(G, 100*milli*darcy, 1);

T = 1*year;
pv = poreVolume(G, rock);
injRate = 10*sum(pv)/T;
nStep = 100;
dT = rampupTimesteps(T, T/nStep);

%%

fluid = initSimpleADIFluid('n'     , [2, 2, 2], ...
                           'rho'   , [1000, 800, 100]*kilogram/meter^3, ...
                           'phases', 'WOG', ...
                           'mu'    , [1, 10, 2]*centi*poise);

model = TwoPhaseOilWaterModel(G, rock, fluid);

W = [];
W = addWell(W, G, rock, 1, 'comp_i', [1, 0], 'type', 'rate', 'val', injRate);
W = addWell(W, G, rock, G.cells.num, 'comp_i', [0, 1], 'type', 'bhp', 'val', 0      );
schedule = simpleSchedule(dT, 'W', W);

%%

state0 = initResSol(G, 100*barsa, [0 1]);
state0.wellSol = initWellSolAD(W, model, state0);

[wsTP, statesTP, reportsTP] = simulateScheduleAD(state0, model, schedule);

%%

fluid = addSolventProperties(fluid, 'n', 2, ...
                                    'rho', 100*kilogram/meter^3, ...
                                    'mixPar', 0, ...
                                    'mu'    , 1*centi*poise);
                                
model = FourPhaseSolventModel(G, rock, fluid);
model.extraStateOutput = true;

W = [];
W = addWell(W, G, rock, 1, 'comp_i', [1, 0, 0, 0], 'type', 'rate', 'val', injRate);
W = addWell(W, G, rock, G.cells.num, 'comp_i', [0, 0, 0, 1], 'type', 'bhp', 'val', 0);
schedule = simpleSchedule(dT, 'W', W);

%%

state0 = initResSol(G, 100*barsa, [0 1 0 0]);
state0.wellSol = initWellSolAD(W, model, state0);

[wsFP, statesFP, reportsFP] = simulateScheduleAD(state0, model, schedule);

%%

close all
figure(1)
plotToolbar(G, statesTP, 'plot1d', true);

figure(2)
plotToolbar(G, statesFP, 'plot1d', true);

%%

p_err = cellfun(@(st1, st2) norm(st1.pressure - st2.pressure)./norm(st1.pressure), statesTP, statesFP);
W_err = cellfun(@(st1, st2) norm(st1.s(:,1) - st2.s(:,1))./norm(st1.s(:,1)), statesTP, statesFP);
O_err = cellfun(@(st1, st2) norm(st1.s(:,2) - st2.s(:,2))./norm(st1.s(:,2)), statesTP, statesFP);
