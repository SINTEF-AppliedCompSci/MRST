mrstModule add dg ad-props ad-core ad-blackoil blackoil-sequential ...
    mrst-gui matlab_bgl vem vemmech vista
mrstVerbose on

%%

setup    = getDGTestCase('simple1d', 'n', 20, 'nkr', 1, 'degree', 0:2, 'useUnstructCubature', true);

%%

% setup.modelFV.transportModel.formulation = 'missingPhase';
[wsFV, stFV, repFV] = simulateScheduleAD(setup.state0, setup.modelFV, setup.schedule);

%%

[wsDG, stDG, repDG] = deal(cell(numel(setup.modelDG),1));

%%

ix = 3;
for i = ix
    setup.modelDG{i}.transportModel.disc.limiter = getLimiter(setup.modelDG{i}.transportModel, 'tvb', 1e-3, 'plot', true);
    [wsDG{i}, stDG{i}, repDG{i}] = simulateScheduleAD(setup.state0, setup.modelDG{i}, setup.schedule);
end

%%

ws = vertcat({wsFV}, wsDG(ix));
datasetNames = horzcat({'FV'}, cellfun(@(ix) ['dG(', num2str(ix-1), ')'], num2cell(ix), 'UniformOutput', false));
plotWellSols(ws, setup.schedule.step.val, 'datasetnames', datasetNames)

%%

close all
t = 50;
figure, hold on
for i = ix
    plotSaturationDG(setup.modelDG{i}.transportModel.disc, stDG{i}{t}, 'plot1d', true, 'n', 500);
end

%%

setupRef = getDGTestCase('simple1d', 'n', 1000, 'nkr', 3, 'degree', 0:5, 'useUnstructCubature', false);
[wsRef, stRef, repRef] = simulateScheduleAD(setupRef.state0, setupRef.modelFV, setupRef.schedule);

%%

ws = vertcat({wsRef}, {wsFV}, wsDG(ix));
datasetNames = horzcat({'Ref', 'FV'}, cellfun(@(ix) ['dG(', num2str(ix-1), ')'], num2cell(ix), 'UniformOutput', false));
plotWellSols(ws, setup.schedule.step.val, 'datasetnames', datasetNames)