mrstModule add dg ad-props ad-core ad-blackoil blackoil-sequential ...
    mrst-gui matlab_bgl vemmech vista
mrstVerbose on

%%

setup   = getDGTestCase('simple1d', 'n', 20);

modelDG = TransportOilWaterDGModel(setup.modelFV);

%%

[wsFV, stFV, repFV] = simulateScheduleAD(setup.state0, setup.modelFV, setup.schedule);

%%

[wsDG, stDG, repDG] = simulateScheduleAD(setup.state0, modelDG, setup.schedule);