mrstModule add dg ad-props ad-core ad-blackoil blackoil-sequential ...
    mrst-gui matlab_bgl vem vemmech vista upr compositional spe10
mrstVerbose on

%%

setup = getDGTestCase('simple1d', 'n', 20, 'nkr', 3); %#ok

%%

setup = getDGTestCase('qfs_wo_2d', 'useMomentFitting', true, 'nkr', 1, 'k', {[0,0], [0,0; 1,0; 0,1; 1,1]}, 'degree', 1); %#ok

%%

setup = getDGTestCase('qfs_wog_3d', 'useMomentFitting', true, 'degree', [0,0,0; 1,1,0; 1,1,1]); %#ok

%%

setup = getDGTestCase('spe10_wo', 'useMomentFitting', true, 'I', 1:30, 'J', 1:110); %#ok

%%

setup = getDGTestCase('spe1', 'subset', false, 'useMomentFitting', true); %#ok

%%

setup = getDGTestCase('qfs_co2_2d', 'n', 3, 'useOverall', false);

%%

sim = @(model,inx) simulateScheduleAD(setup.state0, setup.(model){inx}, setup.schedule);
[wsDG, stDG, repDG] = deal(cell(numel(setup.modelDG),1));

%%

[wsFV, stFV, repFV] = sim('modelFV', 1);

%%

[wsFI, stFI, repFI] = sim('modelFI', 1);

%%

ix = 2;
for i = ix
    [wsDG{i}, stDG{i}, repDG{i}] = sim('modelDG', i);
end

%%

plotWellSols(vertcat(wsDG(ix)), setup.schedule.step.val);

%%

coords = getPlotCoordinates(setup.modelFV{1}.G, 'n', 100);

%%

close all
figure('Position', [0,0,1000,500])

%%

for t = 1:numel(stDG{ix(1)})
    clf
    for i = 1:numel(ix)
        subplot(1, numel(ix), i)
        disc = setup.modelDG{ix(i)}.transportModel.disc;
        plotSaturationDG(disc, stDG{ix(i)}{t}, 'edgecolor', 'none', 'edgealpha', 0.2, 'coords', coords, 'phaseNo', 1);
        axis tight
        zlim([0,1])
        pbaspect([1,1,0.25]);
        view([100,50]);
        camlight
        lighting gouraud
    end
    pause(0.2);
end

%%

ix = 2;
close all
figure('Position', [0,0,1000,500]);
azel = [64,-10];
view(azel);
disc = setup.modelDG{ix(1)}.transportModel.disc;
h = plotSaturationDG(disc, stDG{ix(1)}{1}, 'edgecolor', 'k', 'edgealpha', 0.2, 'coords', coords, 'phaseNo', 3);
for t = 1:numel(stDG{ix(1)})
    for i = 1:numel(ix)
        delete(h);
%         h.Position = [-h.Position(1:2), h.Position(3)];
%         h.Position = -h.Position;
        subplot(1, numel(ix), i)
        h = plotSaturationDG(disc, stDG{ix(i)}{t}, 'edgecolor', 'k', 'edgealpha', 0.2, 'coords', coords, 'phaseNo', 3);
        ax = gca;
        plotWell(setup.modelFV{1}.G, setup.schedule.control(1).W);
%         plotGrid(setup.modelFV{1}.G, 'facec', 'none', 'facealpha', 0.2);
        ax.ZDir = 'reverse';
        box on
        view(azel);
        if t == 1
            hl = camlight;
        end
%         pbaspect([1,1,0.5]);
    end
    pause(0.1);
end