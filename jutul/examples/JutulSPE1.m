%% Example demonstrating SPE1 in Jutul
% The first SPE project comparing black-oil reservoir simulators was
% organized by Odeh (1981) and describes a depletion problem with gas
% injection in a small 10x10x3 reservoir with a producer and an injector
% placed in diagonally opposite corners. The porosity is uniform and equal
% 0.3, whereas the permeability is isotropic with values 500, 50, and 200
% md in the three layers with thickness 20, 30, and 50 ft. The reservoir is
% initially undersaturated with a pressure field that is constant in each
% layer, a uniform mixture of water (Sw = 0.12) and oil (So = 0.88) with no
% initial free gas (Sg = 0.0) and a constant dissolved gas-oil ratio (Rs )
% throughout the model.
%
% This example will set up the problem in MRST, as specified using an
% industry-standard input format that only considers the first 1216 days of
% the original 10-year production horizon. Once set up, we simulate the
% case in both Jutul and MRST.
%
% Odeh, A.S. 1981. Comparison of Solutions to a Three-Dimensional Black-Oil
% Reservoir Simulation Problem. J Pet Technol 33 (1): 13–25. SPE-9723-PA.
% http://dx.doi.org/10.2118/9723-PA 
%

%% Load MRST modules and set up the case
mrstModule add ad-core ad-blackoil deckformat jutul
pth = getDatasetPath('spe1');
deck_path  = fullfile(pth, 'BENCH_SPE1.DATA');
problem = initEclipsePackedProblemAD(deck_path);
schedule = problem.SimulatorSetup.schedule;
model = problem.SimulatorSetup.model;

%% Simulate in Jutul
% The simulation does not start automatically. Instead, the routine
% outputs a command which you must paste into an existing Julia session
% that already has JutulDarcy preloaded. Once this simulation is finished,
% you can hit the return button to continue execution in MATLAB. Notice
% that the first time you run the simulation, it will take a long time
% since Julia has to compile necessary code. Once compiled, however, the
% simulation is fast.
[ws, states] = simulatePackedProblemJutul(problem, 'daemon', false);

%% Simulate in MRST
simulatePackedProblem(problem);
[ws_m, states_m] = getPackedSimulatorOutput(problem);

%% Compare surface gas rates for Jutul and MRST
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
axis off tight, view(50, 50), colorbar('horiz')
title('S_g MRST')

subplot(1, 3, 2)
plotCellData(G, sg)
axis off tight, view(50, 50), colorbar('horiz')
title('S_g Jutul')

subplot(1, 3, 3)
plotCellData(G, sg - sg_m)
axis off tight, view(50, 50), colorbar('horiz')
title('Difference')

%% Compare well responses interactively
n = numel(ws);
T = cumsum(schedule.step.val);
plotWellSols({ws(1:n), ws_m(1:n)}, T(1:n), 'datasetnames', {'Jutul', 'MRST'})