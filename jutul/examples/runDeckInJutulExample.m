mrstModule add ad-core ad-blackoil deckformat ad-props test-suite coarsegrid jutul
if ~exist('name', 'var')
    name = 'spe1';
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
[state0, model, schedule, nonlinear] = initEclipseProblemAD(deck_path, 'ReorderStrategy', reorder);
%% Write case to temporary directory
jpth = writeJutulInput(state0, model, schedule, name);
%% Once simulated, read back as MRST format
[ws, states] = readJutulOutput(jpth);
%%
model.OutputStateFunctions = {'ComponentTotalMass', 'Density', 'PhasePressures', 'RelativePermeability', 'ShrinkageFactors'};
[ws_m, states_m] = simulateScheduleAD(state0, model, schedule, 'NonLinearSolver', getNonLinearSolver(model));
%%
