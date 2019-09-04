mrstModule add dg ad-props ad-core ad-blackoil blackoil-sequential ...
    mrst-gui matlab_bgl vemmech vista
mrstVerbose on

%%

setup = getDGTestCase('simple1d', 'n', 20);

%%

[wsFV, stFV, repFV] = simulateScheduleAD(setup.state0, setup.modelFV, setup.schedule);

%%

[wsDG, stDG, repDG] = deal(cell(numel(setup.modelDG),1));
for dNo = 1:numel(setup.modelDG)
    [wsDG{dNo}, stDG{dNo}, repDG{dNo}] = simulateScheduleAD(setup.state0, setup.modelDG{dNo}, setup.schedule);
end

%%

close all

nc = setup.modelFV.G.cells.num;
x = linspace(1/nc, 1 - 1/nc, nc);

ix = 50;
hold on
plot(x, stFV{ix}.s(:,1));
marker = {'o', 'sq', '^'};
for dNo = 1:numel(setup.modelDG)
    plot(x, stDG{dNo}{ix}.s(:,1), marker{dNo});
end
hold off

%%

close all
dNo = 3;
plotToolbar(setup.modelFV.G, stDG{dNo}, 'plot1d', true);