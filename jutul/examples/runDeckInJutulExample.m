%% Example demonstrating how to run black oil cases from .DATA files in Jutul
mrstModule add ad-core ad-blackoil ad-props deckformat jutul
if ~exist('name', 'var')
    name = 'spe1';
end
if ~exist('use_daemon', 'var')
    use_daemon = false;
end
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
[state0, model, schedule, nls] = initEclipseProblemAD(deck_path, 'ReorderStrategy', reorder);
%% Write case to temporary directory and run
% You can set daemon mode to true if set up, see backgroundJutulExample.
[ws, states] = simulateScheduleJutul(state0, model, schedule, 'daemon', use_daemon);
%% Run in MRST
[ws_m, states_m] = simulateScheduleAD(state0, model, schedule, 'NonLinearSolver', nls);
%% Compare wells
plotWellSols({ws, ws_m}, 'datasetnames', {'Jutul', 'MRST'})