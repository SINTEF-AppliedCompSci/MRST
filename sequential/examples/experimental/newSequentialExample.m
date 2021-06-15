mrstModule add ad-blackoil ad-core mrst-gui ad-props deckformat sequential linearsolvers
if ~exist('basename', 'var')
    basename = 'spe1';
end

if ~exist('useFastTransport', 'var')
    useFastTransport = true;
end

switch basename
    case 'spe9'
        pth = getDatasetPath('spe9');
        fn  = fullfile(pth, 'BENCH_SPE9.DATA');
    case 'spe1'
        pth = getDatasetPath('spe1');
        fn  = fullfile(pth, 'BENCH_SPE1.DATA');
    case 'egg'
        pth = mrstPath('test-datasets');
        fn  = fullfile(pth, 'ad-testdata', 'external', 'TUDelft-EGG', 'BENCH_EGG.DATA');
end

[state0, model, schedule, nls] = initEclipseProblemAD(fn, 'useMex', true, 'rowMajorAD', true);
packer = @(model, name) packSimulationProblem(state0, model, schedule, basename, 'name', name, 'NonLinearSolver', nls);

%% Fully-implicit
problem_fi = packer(model, 'FIM');
simulatePackedProblem(problem_fi, 'restartStep', 1);

%% Sequential
pmodel = PressureModel(model, 'reductionStrategy', 'numerical'); % Numerical strategy
tmodel = TransportModel(model);

psolver = AMGCLSolverAD('tolerance', 1e-4, 'verbose', false);
if useFastTransport && isa(model.AutoDiffBackend, 'DiagonalAutoDiffBackend')
    % Use faster assembly since there's no wells in transport
    tsolver = AMGCLSolverBlockAD('preconditioner', 'relaxation', ...
                                 'relaxation', 'ilu0', 'tolerance', 1e-4, ...
                                 'verbose', false);
    tmodel.parentModel.AutoDiffBackend.deferredAssembly = true;
else
    tsolver = selectLinearSolverAD(model, 'useCPR', false, 'verbose', false);
end
seq = SequentialPressureTransportModel(pmodel, tmodel);
seq.transportNonLinearSolver.LinearSolver = tsolver;
seq.transportNonLinearSolver.maxIterations = 18;  % 25 is too much
seq.transportNonLinearSolver.maxTimestepCuts = 0; % Go back to pressure if transport fails
seq.transportNonLinearSolver.continueOnFailure = true; % Let the outer loop handle failure
seq.transportNonLinearSolver.errorOnFailure = false;

seq.pressureNonLinearSolver.LinearSolver = psolver;

problem_s = packer(seq, 'Seq');
simulatePackedProblem(problem_s, 'restartStep', 1);
%% Sequential fully-implicit (outer loop)
seqfi = seq;
seqfi.stepFunctionIsLinear = false;
seqfi.pressureModel.verbose = false;
seqfi.transportModel.verbose = false;
seqfi.transportNonLinearSolver.verbose = false;
seqfi.pressureNonLinearSolver.verbose = false;

problem_sfi = packer(seqfi, 'SFI');
simulatePackedProblem(problem_sfi, 'restartStep', 1);
%% Plot
problems = {problem_fi, problem_s, problem_sfi};
[ws, states, reports, names, time] = getMultiplePackedSimulatorOutputs(problems,...
                                    'readFromDisk', false, 'readWellSolsFromDisk', true);

plotWellSols(ws, time, 'datasetnames', names);