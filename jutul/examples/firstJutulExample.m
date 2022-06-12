%% First introductory example to Jutul as a MRST accelerator
% JutulDarcy is a reservoir simulator that can be used to accelerate MRST.
% We can set up a case in MRST, and run it in Jutul for increased
% computational performance. Jutul supports black-oil, immiscible and
% compositional models with multisegment wells.
mrstModule add ad-core ad-blackoil spe10 deckformat ad-props test-suite compositional jutul
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
        % Assume some test-suite case
        setup = eval(name);
end
[state0, model, schedule] = deal(setup.state0, setup.model, setup.schedule);
%% Write case to disk
pth = writeJutulInput(state0, model, schedule, name);
disp('Pausing - run the command in Julia and hit any key to continue')
pause()
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