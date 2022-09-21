%% Example demonstrating SPE9 run in Jutul
mrstModule add ad-core ad-blackoil deckformat jutul

pth = getDatasetPath('spe9');
deck_path  = fullfile(pth, 'BENCH_SPE9.DATA');
problem = initEclipsePackedProblemAD(deck_path);

schedule = problem.SimulatorSetup.schedule;
model = problem.SimulatorSetup.model;
%% Simulate in Jutul
[ws, states] = simulatePackedProblemJutul(problem, 'daemon', false);
%% Simulate in MRST
simulatePackedProblem(problem);
[ws_m, states_m] = getPackedSimulatorOutput(problem);
%% Plot states interactively
mrstModule add mrst-gui
G = model.G;
figure; plotToolbar(G, states, 'field', 'pressure')
view(30, 45)
title('Jutul')
figure; plotToolbar(G, states_m, 'field', 'pressure')
view(30, 45)
title('MRST')
%% Plot wells interactively
n = numel(ws);
T = cumsum(schedule.step.val);
plotWellSols({ws(1:n), ws_m(1:n)}, T(1:n), 'datasetnames', {'Jutul', 'MRST'})
%% Compare wells
addir = mrstPath('ad-blackoil');
compare = fullfile(addir, 'examples', 'spe9', 'compare');
smry = readEclipseSummaryUnFmt(fullfile(compare, 'SPE9'));

compd = 1:(size(smry.data, 2));
Tcomp =  smry.get(':+:+:+:+', 'YEARS', compd);
T = convertTo(cumsum(schedule.step.val), year);
C = lines(3);

juplot = @(data) plot(T, data, '--', 'linewidth', 3, 'color', C(1, :));
mrstplot = @(data) plot(T, data, '-', 'linewidth', 1, 'color', C(1, :));
compplot = @(data) plot(Tcomp, data, 'ro', 'linewidth', 2, 'color', C(2, :));

figure;
names = {'PROD13', 'PROD18'};
nn = numel(names);
for i = 1:nn

    name = names{i};

    comp = convertFrom(smry.get(name, 'WBHP', compd), psia)';
    mrst = getWellOutput(ws_m, 'bhp', name);
    jutul = getWellOutput(ws, 'bhp', name);

    subplot(nn, 1, i)
    hold on
    mrstplot(mrst);
    compplot(comp);
    juplot(jutul);
    title(name)
    axis tight
    grid on

    xlabel('Time (years)')
    ylabel('Pressure (Pa)')
end
legend({'MRST', 'ECLIPSE', 'Jutul'})