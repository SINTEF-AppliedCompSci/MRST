mrstModule add ad-core ad-blackoil deckformat ad-props test-suite coarsegrid jutul
name = 'spe1';
reorder = {};
switch name
    case 'olympus_1_reduced'
        pth = fullfile(mrstDataDirectory(), 'olympus', 'ECLIPSE', 'OLYMPUS_1_REDUCED');
        deck_path = fullfile(pth, 'OLYMPUS_1_REDUCED.DATA');
        reorder = 'none';
    case 'spe1'
        pth = getDatasetPath('spe1');
        deck_path  = fullfile(pth, 'BENCH_SPE1.DATA');
    case 'qfs_wo'
        setup = qfs_wo('nkr', 2);
    otherwise
        error('No such case.')
end
[state0, model, schedule, nonlinear] = initEclipseProblemAD(deck_path, 'ReorderStrategy', reorder);
%% Write case to temporary directory
pth = writeJutulInput(state0, model, schedule, name);
%% Once simulated, read back as MRST format
[wells, states] = readJutulOutput(pth);
%%
[ws, states_m] = simulateScheduleAD(state0, model, schedule);
%%
qg = abs(getWellOutput(ws, 'qGs', 'PRODUCER'));
figure(1); clf; hold on
plot(wells.TIME, qg, 'linewidth', 2)
plot(wells.TIME, abs(wells.PRODUCER.GRAT), '.', 'MarkerSize', 18)
legend('MRST', 'Jutul')
%%
mrstModule add mrst-gui
G = model.G;
figure; plotToolbar(G, states)
title('Jutul')
figure; plotToolbar(G, states_m)
title('MRST')

