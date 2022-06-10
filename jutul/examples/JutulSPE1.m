%% Example demonstrating SPE1 in Jutul
pth = getDatasetPath('spe1');
deck_path  = fullfile(pth, 'BENCH_SPE1.DATA');
problem = initEclipsePackedProblemAD(deck_path);
schedule = problem.SimulatorSetup.schedule;
model = problem.SimulatorSetup.model;
%% Simulate in Jutul
[ws, states] = simulatePackedProblemJutul(problem);
%% Simulate in MRST
simulatePackedProblem(problem);
[ws_m, states_m] = getPackedSimulatorOutput(problem);
%% Compare a 
t = cumsum(schedule.step.val);
qg = abs(getWellOutput(ws, 'qGs', 'PRODUCER'));
qg_m = abs(getWellOutput(ws_m, 'qGs', 'PRODUCER'));

figure(1); clf; hold on
plot(t, qg_m, 'linewidth', 2)
plot(t, qg, '.', 'MarkerSize', 18)
legend('MRST', 'Jutul')
%% Plot the final gas saturation
sg = states{end}.s(:, 3);
sg_m = states_m{end}.s(:, 3);
mrstModule add mrst-gui
G = model.G;
f = figure();
fp = get(f, 'DefaultFigurePosition');
set(f, 'Position', fp.*[1, 1, 2, 1]);
subplot(1, 3, 1)
plotCellData(G, sg_m)
axis off tight
view(50, 50)
title('S_g MRST')
colorbar('horiz')

subplot(1, 3, 2)
plotCellData(G, sg)
axis off tight
view(50, 50)
title('S_g Jutul')
colorbar('horiz')

subplot(1, 3, 3)
plotCellData(G, sg - sg_m)
axis off tight
view(50, 50)
colorbar('horiz')
title('Difference')
%% Compare wells
n = numel(ws);
T = cumsum(schedule.step.val);
plotWellSols({ws(1:n), ws_m(1:n)}, T(1:n), 'datasetnames', {'Jutul', 'MRST'})