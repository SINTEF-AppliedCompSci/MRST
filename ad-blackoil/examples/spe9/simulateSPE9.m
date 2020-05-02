%% Minimal example that runs SPE9
mrstModule add ad-blackoil ad-core mrst-gui ad-props deckformat
pth = getDatasetPath('spe9');
fn  = fullfile(pth, 'BENCH_SPE9.DATA');
problem = initEclipsePackedProblemAD(fn, 'useMex', true, 'rowMajorAD', true);
%% Simulate
simulatePackedProblem(problem, 'restartStep', 1);
%% Plot
plotPackedProblem(problem);
