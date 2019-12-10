mrstModule add dg ad-props ad-core ad-blackoil blackoil-sequential ...
    mrst-gui matlab_bgl vem vemmech vista upr
mrstVerbose on

%%

setup = getDGTestCase('qfs_wog_3d', 'useUnstructCubature', true, 'n', 3, 'pebi', true);

%%

sim = @(model,inx) simulateScheduleAD(setup.state0, setup.(model){inx}, setup.schedule);

%%

[wsFV, stFV, repFV] = sim('modelFV', 1);

%%

[wsFI, stFI, repFI] = sim('modelFI', 1);

%%

[wsDG, stDG, repDG] = deal(cell(numel(setup.modelDG),1));

%%

[wsDG{1}, stDG{1}, repDG{1}] = sim('modelDG', 1);

%%

[wsDG{2}, stDG{2}, repDG{2}] = sim('modelDG', 2);

%%
    
plotWellSols({wsFV, wsDG{2}});