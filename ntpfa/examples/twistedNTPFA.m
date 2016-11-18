gravity off
mrstModule add ntpfa mimetic incomp

dims = [15, 15];
pdims = [1, 1];

G = cartGrid(dims, pdims);
G = twister(G, 0.1);
G = computeGeometry(G);
fluid = initSimpleADIFluid();

rock = makeRock(G, 1, 1);
model = PressureOilWaterModelNTPFA(G, rock, fluid);
model.verbose = true;

[bc, src, W] = deal([]);
bc = pside(bc, G, 'xmin', 1, 'sat', [1, 0]);
bc = pside(bc, G, 'xmax', 0, 'sat', [1, 0]);
% bc = [];
% src = addSource(src, 1, 1, 'sat', [1, 0]);
% src = addSource(src, G.cells.num, -1, 'sat', [1, 0]);

%%
% NTPFA
state0 = initResSol(G, 0, [1, 0]);
dT = 1;

schedule = simpleSchedule(dT, 'W', W, 'bc', bc, 'src', src);
[~, states] = simulateScheduleAD(state0, model, schedule);
ntpfa = states{1};
% MIMETIC
S = computeMimeticIP(G, rock);
fluid2 = initSimpleFluid('mu', [1, 1], 'rho', [1, 1], 'n', [1, 1]);
mimetic = incompMimetic(state0, G, S, fluid2, 'bc', bc, 'src', src);

% TPFA
T = computeTrans(G, rock);
fluid2 = initSimpleFluid('mu', [1, 1], 'rho', [1, 1], 'n', [1, 1]);
tpfa = incompTPFA(state0, G, T, fluid2, 'bc', bc, 'src', src);
%%
figure(1), clf
subplot(1, 3, 1);
plotCellData(G, ntpfa.pressure)
title('N-TPFA')

subplot(1, 3, 2);
plotCellData(G, mimetic.pressure)
title('Mimetic')

subplot(1, 3, 3);
plotCellData(G, tpfa.pressure)
title('TPFA')


figure; hold on
x = G.cells.centroids(:, 1);
plot(x, ntpfa.pressure, '.')
plot(x, mimetic.pressure, '.')
plot(x, tpfa.pressure, '.')
legend('NTPFA', 'Mimetic', 'TPFA')