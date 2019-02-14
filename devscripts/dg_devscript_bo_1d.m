mrstModule add dg ad-core ad-props ad-blackoil blackoil-sequential weno
mrstVerbose on

%%

n     = 10;
G     = computeGeometry(cartGrid([n, 1], [1000, 100]));
G     = computeCellDimensions2(G);
rock  = makeRock(G, 100*milli*darcy, 0.4);
fluid = initSimpleADIFluid();

modelFI  = ThreePhaseBlackOilModel(G, rock, fluid);
modelSI = getSequentialModelFromFI(modelFI);
modelSI.transportModel.conserveWater = true;


transportModel = TransportBlackOilModelDG(G, rock, fluid, 'degree', 1);
transportModel.conserveWater = true;

modelDG  = SequentialPressureTransportModelDG(modelSI.pressureModel, transportModel);
modelDG.pressureModel.extraWellSolOutput = true;

modelFIWO = TwoPhaseOilWaterModel(G, rock, fluid);
modelSIWO = getSequentialModelFromFI(modelFIWO);
transportModelWO = TransportOilWaterModelDG(G, rock, fluid, 'degree', 1);
transportModelWO.conserveWater = true;
transportModelWO.conserveOil   = true;
modelDGWO = SequentialPressureTransportModelDG(modelSIWO.pressureModel, transportModelWO);

time = 1*year;
rate = 2*sum(poreVolume(G, rock))/time;
bhp = 50*barsa;

W = [];
W = addWell(W, G, rock, 1, 'type', 'rate', 'val', rate, 'name', 'I');
W = addWell(W, G, rock, G.cells.num, 'type', 'bhp' , 'val', bhp , 'name', 'P');

WWO = [];
WWO = addWell(WWO, G, rock, 1, 'type', 'rate', 'val', rate, 'name', 'I', 'comp_i', [1,0]);
WWO = addWell(WWO, G, rock, G.cells.num, 'type', 'bhp' , 'val', bhp , 'name', 'P', 'comp_i', [1,0]);

sO       = 1;
state0   = initResSol(G, bhp, [0, sO, 1-sO]);
state0WO = initResSol(G, bhp, [0, sO]);
% state0 = assignDofFromState(modelDG.transportModel.disc, state0);

dt       = 7*day;
dtvec    = rampupTimesteps(time, dt, 0);
schedule = simpleSchedule(dtvec, 'W', W);
scheduleWO = simpleSchedule(dtvec, 'W', WWO);

%%

[wsDG_bo, stDG_bo, repDG_bo] = simulateScheduleAD(state0, modelDG, schedule);

%%

[wsDG_ow, stDG_ow, repDG_ow] = simulateScheduleAD(state0WO, modelDGWO, scheduleWO);

%%

[wsFV, stFV, repFV] = simulateScheduleAD(state0, modelSI, schedule);

%%

close all
plotWellSols({wsFV, wsDG_bo, wsDG_ow}, schedule.step.val)

%%

close all
plotToolbar(G, stDG)
