mrstModule add dg vem vemmech ad-props ad-core ad-blackoil blackoil-sequential

%%

n = 3;
G = computeGeometry(cartGrid([n,1], [100,10]));
G.nodes.coords = G.nodes.coords;
G = computeVEMGeometry(G);

rock = makeRock(G, 100*milli*darcy, 0.4);
fluid = initSimpleADIFluid('phases', 'WO', ...
                           'rho', [1000, 800]*kilogram/meter^3, ...
                           'mu', [0.3, 1]*centi*poise);
                       
modelfi = TwoPhaseOilWaterModel(G, rock, fluid);
model = getSequentialModelFromFI(modelfi);
model.transportModel = TransportOilWaterModelDG(G, rock, fluid, 'degree', 0);
                       
%%

time = 2*year;
rate = 1*sum(poreVolume(G, rock))/time;
W = [];
W = addWell(W, G, rock, 1, 'type', 'rate', 'val', rate, 'comp_i', [1,0]);
W = addWell(W, G, rock, G.cells.num, 'type', 'bhp', 'val', 50*barsa, 'comp_i', [1,0]);


dt = 30*day;
dtvec = rampupTimesteps(time, dt, 0);

schedule = simpleSchedule(dtvec, 'W', W);

%%

state0      = initResSol(G, 100*barsa, [0,1]);
nDof        = model.transportModel.basis.nDof;
state0.sdof = zeros(G.cells.num*nDof, 2);
state0.sdof(1:nDof:G.cells.num*nDof,2) = 1;

[ws, state, rep] = simulateScheduleAD(state0, model, schedule);

%%

basis = dgBasis(model.transportModel.degree, model.transportModel.G.griddim, 'poly');
x = zeros(G.cells.num,2);

sW = cellfun(@(s) getSatFromDof(x, (1:G.cells.num)', s.sdof(:,1), model.transportModel), state, 'unif', false);

plotToolbar(G, state);