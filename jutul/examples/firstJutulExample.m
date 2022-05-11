mrstModule add ad-core ad-blackoil spe10 deckformat ad-props test-suite coarsegrid compositional jutul
close all; clear
if ~exist('name', 'var')
    name = 'qfs_wo';
end
switch name
    case 'qfs_wo'
        setup = qfs_wo();
    case 'buckley_leverett_wo'
        setup = buckley_leverett_wo();
    case 'spe10_layer'
        setup = spe10_wo('layers', 1);
    case 'fractures_compositional'
        setup = fractures_compositional();
    otherwise
        setup = eval(name);
end
[state0, model, schedule] = deal(setup.state0, setup.model, setup.schedule);
%% Write case to disk
pth = writeJutulInput(state0, model, schedule, name);
%% Once simulated, read back as MRST format
[ws, states] = readJutulOutput(pth);
%% Simulate MRST for comparison purposes
nls = getNonLinearSolver(model);
[ws_m, states_m] = simulateScheduleAD(state0, model, schedule, 'NonLinearSolver', nls);
%% Compare the results
mrstModule add mrst-gui
G = model.G;
figure;
plotToolbar(G, states_m)
title('MRST')
figure;
plotToolbar(G, states)
title('Jutul')
figure;
plotToolbar(G, applyFunction(@(x, y) compareStructs(x, y), states, states_m))
title('Difference')
%% Plot and compare wells
% Jutul uses multisegment wells by default which gives small differences in
% well curves
plotWellSols({ws, ws_m}, cumsum(schedule.step.val), 'datasetnames', {'Jutul', 'MRST'})