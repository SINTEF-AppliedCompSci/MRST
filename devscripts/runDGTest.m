mrstModule add dg ad-props ad-core ad-blackoil blackoil-sequential ...
    mrst-gui matlab_bgl vem vemmech vista upr
mrstVerbose on

%%

setup = getDGTestCase('simple1d', ...
                      'n'                  , 20        , ...
                      'nkr'                , 2         , ...
                      'degree'             , [0,0; 1,0; 2,0; 3,0], ...
                      'basis'              , 'legendre', ...
                      'useUnstructCubature', false     ); %#ok

%%

setup = getDGTestCase('qfs_wo_2d', 'useUnstructCubature', true); %#ok

%%

setup = getDGTestCase('qfs_wog_3d', 'useUnstructCubature', true, 'degree', [1,1,0]);

%%

sim = @(model,inx) simulateScheduleAD(setup.state0, setup.(model){inx}, setup.schedule);
[wsDG, stDG, repDG] = deal(cell(numel(setup.modelDG),1));

%%

[wsFV, stFV, repFV] = sim('modelFV', 1);

%%

ix = 4;
for i = ix
    [wsDG{i}, stDG{i}, repDG{i}] = sim('modelDG', i);
end

%%

plotWellSols(vertcat({wsFV}, wsDG(ix)));

%%

ix = 1:2;

close all
for t = 1:numel(stDG{ix(1)})
    for i = 1:numel(ix)
        subplot(1, numel(ix), i)
        disc = setup.modelDG{ix(i)}.transportModel.disc;
        plotSaturationDG(disc, stDG{ix(i)}{t}, 'edgecolor', 'k', 'edgealpha', 0.2);
        view(3);
        zlim([0,1])
        pbaspect([1,1,0.5]);
        view([100,50]);
    end
    pause(0.1);
end