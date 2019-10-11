n = 20;
setup     = getDGTestCase('qfs_wo_2d', 'useMomentFitting', true, 'nkr', 1, 'n', n);
setup_old = getDGTestCase('qfs_wo_2d', 'useMomentFitting', true, 'nkr', 1, 'n', n, 'useGenericFV', false);

%%

sim = @(setup, model, inx) simulateScheduleAD(setup.state0, setup.(model){inx}, setup.schedule);
[wsDG, stDG, repDG] = deal(cell(numel(setup.modelDG),1));

%%

[wsFV, stFV, repFV] = sim(setup, 'modelFV', 1);

%%

[wsFV_old, stFV_old, repFV_old] = sim(setup_old, 'modelFV', 1);

%%

plotWellSols({wsFV, wsFV_old})

%%

close all
sd = cellfun(@(s1, s2) compareStates(s1,s2), stFV, stFV_old, 'UniformOutput', false);
plotToolbar(setup.modelFV{1}.G, sd);
axis equal tight

%%

setup     = getDGTestCase('spe1');
setup_old = getDGTestCase('spe1', 'useGenericFV', false);

%%

sim = @(setup, model, inx) simulateScheduleAD(setup.state0, setup.(model){inx}, setup.schedule);
[wsDG, stDG, repDG] = deal(cell(numel(setup.modelDG),1));

%%

[wsFV, stFV, repFV] = sim(setup, 'modelFV', 1);

%%

[wsFV_old, stFV_old, repFV_old] = sim(setup_old, 'modelFV', 1);

%%

[wsFI, stFI, repFI] = sim(setup, 'modelFI', 1);

%%

plotWellSols({wsFI, wsFV, wsFV_old})

%%

close all
sd = cellfun(@(s1, s2) compareStates(s1,s2), stFV, stFV_old, 'UniformOutput', false);
plotToolbar(setup.modelFV{1}.G, sd);
axis equal tight

%%

plotToolbar(setup.modelFV{1}.G, stFV);

%%

setup     = getDGTestCase('qfs_co2_2d', 'n', 3, 'useGeneric', true, 'useOverall', false);
setup_old = getDGTestCase('qfs_co2_2d', 'n', 3, 'useGeneric', false, 'useOverall', false);

%%

sim = @(setup, model, inx) simulateScheduleAD(setup.state0, setup.(model){inx}, setup.schedule);
[wsDG, stDG, repDG] = deal(cell(numel(setup.modelDG),1));

%%

[wsFI, stFI, repFI] = sim(setup, 'modelFI', 1);

%%

[wsFI_old, stFI_old, repFI_old] = sim(setup_old, 'modelFI', 1);

%%

plotWellSols({wsFI_old, wsFI})
