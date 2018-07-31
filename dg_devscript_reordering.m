mrstModule add dg vem vemmech ad-props ad-core ad-blackoil blackoil-sequential gasinjection reorder matlab_bgl

%%

n = 3;
l = 1000;
G = computeGeometry(cartGrid([n,n], [l,l]*meter));
G = computeVEMGeometry(G);
G = computeCellDimensions(G);

rock = makeRock(G, 100*milli*darcy, 1);
fluid = initSimpleADIFluid('phases', 'WO'                   , ...
                           'rho'   , [1, 1]*kilogram/meter^3, ...
                           'mu'    , [0.5, 0.5]*centi*poise     , ...
                           'n'     , [1, 1]                 );

modelfi = TwoPhaseOilWaterModel(G, rock, fluid);
modelFV = getSequentialModelFromFI(modelfi);
modelDG = modelFV;

%%

degree = 1;
disc   = DGDiscretization(modelDG.transportModel, 2, 'degree', degree, ...
                         'basis', 'legendre', 'useUnstructCubature', true);
modelDG.transportModel = TransportOilWaterModelDG(G, rock, fluid, 'disc', disc);    

%%

time = 2*year;
rate = 1*sum(poreVolume(G, rock))/time;
W = [];
% W = addWell(W, G, rock, 1          , 'type', 'rate', 'val', rate    , 'comp_i', [1,0]);
W = addWell(W, G, rock, G.cells.num, 'type', 'bhp' , 'val', 50*barsa, 'comp_i', [1,0]);

bc = [];
bc = fluxside(bc, G, 'left', rate, 'sat', [1,0]);
bc = makeDGBC(disc, bc);

dt    = 30*day;
dtvec = rampupTimesteps(time, dt, 0);
schedule = simpleSchedule(dtvec, 'W', W, 'bc', bc);

sW     = 0.0;
state0 = initResSol(G, 100*barsa, [sW,1-sW]);
state0.cells = (1:G.cells.num)';
state0 = assignDofFromState(modelDG.transportModel.disc, state0);

%%

[wsDG, statesDG, rep] = simulateScheduleAD(state0, modelDG, schedule);

%%

[modelDG.transportModel.extraStateOutput, modelDG.pressureModel.extraStateOutput] = deal(true);
modelDGreorder = modelDG;
modelDGreorder.pressureModel.extraStateOutput = true;

modelDGreorder.transportModel = ReorderingModelDG(modelDGreorder.transportModel);

modelDGreorder.transportModel.chunkSize = 1;
modelDGreorder.transportModel.parent.extraStateOutput = true;


[~, states_reorder, rep_reorder] = simulateScheduleAD(state0, modelDGreorder, schedule);

