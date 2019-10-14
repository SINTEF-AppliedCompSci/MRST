mrstModule add dg ad-props ad-core ad-blackoil blackoil-sequential ...
    mrst-gui matlab_bgl vem vemmech vista
mrstVerbose on

%%

setup = getDGTestCase('spe1');

%%

ix = 1:7;
schdl = setup.schedule;
schdl.step.val = schdl.step.val(ix);
schdl.step.control= schdl.step.control(ix);
sim = @(model, inx) simulateScheduleAD(setup.state0, setup.(model){inx}, schdl);

%%

[wsFV, stFV, repFV] = sim('modelFV', 1);

%%

[wsDG{1}, stDG{1}, repDG{1}] = sim('modelDG', 1);

%%

[wsDG{2}, stDG{2}, repDG{2}] = sim('modelDG', 2);

%%

[wsFI, stFI, repFI] = sim('modelFI', 1);

%%

plotWellSols({wsFI, wsDG{:}})