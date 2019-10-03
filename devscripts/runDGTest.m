mrstModule add dg ad-props ad-core ad-blackoil blackoil-sequential ...
    mrst-gui matlab_bgl vem vemmech vista upr
mrstVerbose on

%%

setup = getDGTestCase('simple1d', 'n', 20, 'nkr', 2 ); %#ok

%%

setup = getDGTestCase('qfs_wo_2d', 'useMomentFitting', true); %#ok

%%

setup = getDGTestCase('qfs_wog_3d', 'useMomentFitting', true, 'degree', [1,1,1]);

%%

sim = @(model,inx) simulateScheduleAD(setup.state0, setup.(model){inx}, setup.schedule);
[wsDG, stDG, repDG] = deal(cell(numel(setup.modelDG),1));

%%

[wsFV, stFV, repFV] = sim('modelFV', 1);

%%

ix = 1:2;
for i = ix
    [wsDG{i}, stDG{i}, repDG{i}] = sim('modelDG', i);
end

%%

plotWellSols(vertcat({wsFV}, wsDG(ix)));

%%

coords = getPlotCoordinates(setup.modelFV{1}.G, 'n', 100);

%%

close all

for t = 1:numel(stDG{ix(1)})
    clf
    for i = 1:numel(ix)
        subplot(1, numel(ix), i)
        disc = setup.modelDG{ix(i)}.transportModel.disc;
        plotSaturationDG(disc, stDG{ix(i)}{t}, 'edgecolor', 'none', 'edgealpha', 0.2, 'coords', coords, 'phaseNo', 1);
        axis tight
        zlim([0,1])
        pbaspect([1,1,0.5]);
        view([100,50]);
        camlight
        lighting gouraud
    end
    pause(0.1);
end

%%

ix = 1;
close all
azel = [13,25];
for t = 1:numel(stDG{ix(1)})
    for i = 1:numel(ix)
        clf;
        view(azel);
        camlight;
        h = camlight;
        h.Position = [-h.Position(1:2), h.Position(3)];
        subplot(1, numel(ix), i)
        disc = setup.modelDG{ix(i)}.transportModel.disc;
        plotSaturationDG(disc, stDG{ix(i)}{t}, 'edgecolor', 'k', 'edgealpha', 0.2, 'coords', coords, 'phaseNo', 3);
        ax = gca;
        ax.ZDir = 'reverse';
%         pbaspect([1,1,0.5]);
    end
    pause(0.1);
end