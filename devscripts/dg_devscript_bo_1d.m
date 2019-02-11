mrstModule add dg ad-core ad-props ad-blackoil blackoil-sequential weno
mrstVerbose on

%%

n     = 10;
G     = computeGeometry(cartGrid([n, 1], [1000, 100]));
G     = computeCellDimensions2(G);
rock  = makeRock(G, 100*milli*darcy, 0.4);
fluid = initSimpleADIFluid();
fluid.pcOW = @(s) (1-s)*barsa;

modelFI  = ThreePhaseBlackOilModel(G, rock, fluid);
modelSI = getSequentialModelFromFI(modelFI);
modelSI.transportModel.conserveWater = true;


transportModel = TransportBlackOilModelDG(G, rock, fluid, 'degree', 1);
transportModel.conserveWater = true;
modelDG  = SequentialPressureTransportModelDG(modelSI.pressureModel, transportModel);
modelDG.pressureModel.extraWellSolOutput = true;

time = 1*year;
rate = sum(poreVolume(G, rock))/time;
bhp = 50*barsa;

W = [];
W = addWell(W, G, rock, 1, 'type', 'rate', 'val', rate, 'name', 'I');
W = addWell(W, G, rock, G.cells.num, 'type', 'bhp' , 'val', bhp , 'name', 'P');

sO     = 0.5;
state0 = initResSol(G, bhp, [0, sO, 1-sO]);
% state0 = assignDofFromState(modelDG.transportModel.disc, state0);

dt       = 7*day;
dtvec    = rampupTimesteps(time, dt, 0);
schedule = simpleSchedule(dtvec, 'W', W);

%%

[wsDG, stDG, repDG] = simulateScheduleAD(state0, modelDG, schedule);

%%

[wsFV, stFV, repFV] = simulateScheduleAD(state0, modelSI, schedule);

%%

plotWellSols({wsFV, wsDG}, schedule.step.val)

%%

close all
plotToolbar(G, stDG)

