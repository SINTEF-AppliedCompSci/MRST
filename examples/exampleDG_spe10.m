mrstModule add spe10

%%

[state0, model, schedule]  = setupSPE10_AD('layers', 10);
G = model.G;
fluid = model.fluid;
rock = model.rock;

modelFV = getSequentialModelFromFI(model);
modelDG = modelFV;
modelDG.transportModel = TransportOilWaterModelDG(G, rock, fluid);    
disc    = DGDiscretization(modelDG.transportModel, 2, 'degree', 1, 'basis', 'legendre');

modelDG.transportModel.disc = disc;

state0 = disc.assignDofFromState(state0);
[wsDG, statesDG, rep] = simulateScheduleAD(state0, modelDG, schedule);

