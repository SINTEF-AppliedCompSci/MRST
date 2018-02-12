mrstModule add dg vem vemmech ad-props ad-core ad-blackoil blackoil-sequential

%%

n = 2;
G = computeGeometry(cartGrid([n,n], [2,2]));
G.nodes.coords = G.nodes.coords;
G = computeVEMGeometry(G);

rock = makeRock(G, 100*milli*darcy, 0.4);
fluid = initSimpleADIFluid('phases', 'WO', ...
                           'rho', [1000, 800]*kilogram/meter^3, ...
                           'mu', [0.3, 1]*centi*poise);
                       
modelfi = TwoPhaseOilWaterModel(G, rock, fluid);
model = getSequentialModelFromFI(modelfi);
model.transportModel = TransportOilWaterModelDG(G, rock, fluid);
                       
%%

time = 2*year;
rate = 0.2*sum(poreVolume(G, rock))/time;
W = [];
W = addWell(W, G, rock, 1, 'type', 'rate', 'val', rate, 'comp_i', [1,0]);
W = addWell(W, G, rock, G.cells.num, 'type', 'bhp', 'val', 50*barsa, 'comp_i', [1,0]);


dt = 30*day;
dtvec = rampupTimesteps(time, dt, 0);

schedule = simpleSchedule(dtvec, 'W', W);

%%

state0 = initResSol(G, 100*barsa, [1,0]);
[k, nDof] = dgBasis(model.transportModel.degree, G.griddim);
state0.sdof = zeros(G.cells.num*nDof, 1);
state0.sdof(1:nDof:G.cells.num*nDof) = 1;
[ws, state, rep] = simulateScheduleAD(state0, model, schedule);