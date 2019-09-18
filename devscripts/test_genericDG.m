mrstModule add dg ad-props ad-core ad-blackoil blackoil-sequential ...
    mrst-gui matlab_bgl vemmech vista
mrstVerbose on

%%

setup   = getDGTestCase('simple1d', 'n', 20);

%%

% setup.modelFV.transportModel.formulation = 'missingPhase';
[wsFV, stFV, repFV] = simulateScheduleAD(setup.state0, setup.modelFV, setup.schedule);

%%

[wsDG, stDG, repDG] = deal(cell(numel(setup.modelDG),1));

%%

ix = 1:3;
for i = ix
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