mrstModule add dg ad-props ad-core ad-blackoil blackoil-sequential ...
    mrst-gui matlab_bgl vem vemmech vista upr
mrstVerbose on

%%

setup111 = getDGTestCase('qfs_wog_3d', 'useUnstructCubature', true, 'degree', [1,1,1], 'n', 10);
setup110 = getDGTestCase('qfs_wog_3d', 'useUnstructCubature', true, 'degree', [1,1,0], 'n', 10);

%%

[wsFV, stFV, repFV] = simulateScheduleAD(setup111.state0, setup111.modelFV{1}, setup111.schedule);

%%

[wsDG111, stDG111, rep111] = simulateScheduleAD(setup111.state0, setup111.modelDG{1}, setup111.schedule);

%%

[wsDG110, stDG110, rep110] = simulateScheduleAD(setup110.state0, setup110.modelDG{1}, setup110.schedule);

%%

plotWellSols({wsDG111, wsDG111})