mrstModule add dg ad-props ad-core ad-blackoil blackoil-sequential ...
    mrst-gui matlab_bgl vem vemmech vista
mrstVerbose on

%%

setup    = getDGTestCase('qfs_wo_2d', 'n', 10, 'nkr', 1);
setup.modelDG{2}.transportModel.formulation = 'missingPhase';
sim = @(model, inx) simulateScheduleAD(setup.state0, setup.(model){inx}, setup.schedule);

%%

[wsFV, stFV, repFV] = sim('modelFV', 1);

%%

[wsDG, stDG, repDG] = sim('modelDG', 3);

%%

[wsFI, stFI, repFI] = sim('modelFI', 1);

%%

plotWellSols({wsFI, wsFV, wsDG});

%%

ix = 3;
disc = setup.modelDG{ix}.transportModel.disc;

close all
for t = 1:numel(stDG)
    plotSaturationDG(disc, stDG{t}, 'edgecolor', 'k', 'edgealpha', 0.2);
    view(3);
    zlim([0,1])
    view([100,50]);
    pause(0.1);
end

%%

G = setup.modelFI{1}.G;
plotToolbar(G, stDG)