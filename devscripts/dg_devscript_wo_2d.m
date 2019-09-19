mrstModule add dg ad-props ad-core ad-blackoil blackoil-sequential ...
    mrst-gui matlab_bgl vem vemmech vista
mrstVerbose on

%%

setup    = getDGTestCase('qfs_wo_2d', 'n', 3);
sim = @(model, inx) simulateScheduleAD(setup.state0, setup.(model){inx}, setup.schedule);

%%

[wsFV, stFV, repFV] = sim('modelFV', 1);

%%

[wsDG, stDG, repDG] = sim('modelDG', 1);

%%

[wsFI, stFI, repFI] = sim('modelFI', 1);

%%

plotWellSols({wsFI, wsFV, wsDG});