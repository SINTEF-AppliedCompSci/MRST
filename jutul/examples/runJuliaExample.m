mrstModule add ad-core ad-blackoil spe10 deckformat ad-props test-suite coarsegrid compositional jutul
close all; clear
name = 'qfs_wo';
% name = 'spe10_layer'
% name = 'fractures_compositional';
% name = 'spe10_wo'
% name = '1d_validation'
setup = [];
switch name
    case 'qfs_wo'
        setup = qfs_wo('nkr', 2, 'ncells', 2, 'pvi', 1);
    case 'buckley_leverett_wo'
        setup = buckley_leverett_wo();
    case 'spe10_layer'
        setup = spe10_wo('layers', 1);
    case '1d_validation'
        [state0, model, schedule, ref] = setupSimpleCompositionalExample(false);
    otherwise
        setup = eval(name);
end
if ~isempty(setup)
    % Was test case of some sorts
    [state0, model, schedule] = deal(setup.state0, setup.model, setup.schedule);
end
%%
pth = writeJutulInput(state0, model, schedule, name);
%% Once simulated, read back as MRST format
[wells, states] = readJutulOutput(pth);
%%
nls = getNonLinearSolver(model, 'TimestepStrategy', 'none');
[ws, states_m] = simulateScheduleAD(state0, model, schedule, 'NonLinearSolver', nls);
%%
close all
mrstModule add mrst-gui
G = model.G;
%%
figure;
plotToolbar(G, states_m)
title('MRST')
figure;
plotToolbar(G, states)
title('Jutul')
figure;
plotToolbar(G, applyFunction(@(x, y) compareStructs(x, y), states, states_m))
title('Difference')