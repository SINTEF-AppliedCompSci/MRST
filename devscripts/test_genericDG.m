mrstModule add dg ad-props ad-core ad-blackoil blackoil-sequential ...
    mrst-gui matlab_bgl vemmech vista
mrstVerbose on

%%

setup   = getDGTestCase('simple1d', 'n', 20);

%%

[wsFV, stFV, repFV] = simulateScheduleAD(setup.state0, setup.modelFV, setup.schedule);

%%

for dNo = 2:numel(setup.modelDG)
    [wsDG, stDG, repDG] = simulateScheduleAD(setup.state0, setup.modelDG{dNo}, setup.schedule);
end