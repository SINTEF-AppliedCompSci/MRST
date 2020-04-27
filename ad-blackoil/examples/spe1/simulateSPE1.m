%% Minimal example that runs SPE1
mrstModule add ad-blackoil ad-core mrst-gui ad-props deckformat
pth = getDatasetPath('spe1');
fn  = fullfile(pth, 'BENCH_SPE1.DATA');
problem = initEclipsePackedProblemAD(fn);
%% Simulate
simulatePackedProblem(problem, 'restartStep', 1);
%% Plot
plotPackedProblem(problem);
