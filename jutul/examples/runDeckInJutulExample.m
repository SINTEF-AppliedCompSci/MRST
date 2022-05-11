mrstModule add ad-core ad-blackoil deckformat ad-props test-suite coarsegrid jutul
name = 'spe1';
reorder = {};
switch name
    case 'spe1'
        pth = getDatasetPath('spe1');
        deck_path  = fullfile(pth, 'BENCH_SPE1.DATA');
    case 'spe9'
        pth = getDatasetPath('spe9');
        deck_path  = fullfile(pth, 'BENCH_SPE9.DATA');
    case 'egg'
        deck_path = getDeckEGG();
    otherwise
        error('No such case.')
end
[state0, model, schedule, nonlinear] = initEclipseProblemAD(deck_path, 'ReorderStrategy', reorder);
%% Write case to temporary directory
pth = writeJutulInput(state0, model, schedule, name);
%% Once simulated, read back as MRST format
[ws, states] = readJutulOutput(pth);
%%
[ws_m, states_m] = simulateScheduleAD(state0, model, schedule);
%%
t = cumsum(schedule.step.val);
qg = abs(getWellOutput(ws, 'qGs', 'PRODUCER'));
qg_m = abs(getWellOutput(ws_m, 'qGs', 'PRODUCER'));

figure(1); clf; hold on
plot(t, qg_m, 'linewidth', 2)
plot(t, qg, '.', 'MarkerSize', 18)
legend('MRST', 'Jutul')
%%
mrstModule add mrst-gui
G = model.G;
figure; plotToolbar(G, states)
title('Jutul')
figure; plotToolbar(G, states_m)
title('MRST')

