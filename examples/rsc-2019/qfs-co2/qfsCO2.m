mrstModule add spe10 vem vemmech dg ad-core ad-props ad-blackoil ...
    blackoil-sequential incomp msrsb matlab_bgl coarsegrid mrst-gui ...
    reorder weno vista compositional

mrstVerbose on

%% Common stuff

name = 'qfs_co2_small';
baseName = 'qfs-CO2';
dataDir  = fullfile('/media/strene/806AB4786AB46C92/mrst-dg/rsc-2019', baseName);
[state0, model, schedule] = getReorderingCase(name);

pack = @(model, name, nls, desc) ...
    packSimulationProblem(state0, model, schedule, baseName, ...
                          'Name'           , name          , ...
                          'NonLinearSolver', nls           , ...
                          'Directory'      , dataDir       , ...
                          'Description'    , desc          );


%% Fully-implicit

nls = NonLinearSolver();
nls.LinearSolver = selectLinearSolverAD(model);
model.AutoDiffBackend = DiagonalAutoDiffBackend();

fim = pack(model, 'fim', nls, 'Fully-implicit');

%% Sequential-implicit

modelseq = getSequential_local(model);
nls = NonLinearSolver();
seq = pack(modelseq, 'seq', nls, 'Sequential');

%% Sequential-reordering

nls = NonLinearSolver();

modelreorder = getSequential_local(model);
modelreorder.pressureModel.extraStateOutput = true;
modelreorder.transportModel = ReorderingModel(modelreorder.transportModel);
modelreorder.transportModel.parent.extraStateOutput = true;

modelreorder.transportModel.chunkSize = 1;
modelreorder.transportModel.buffer = 0;
modelreorder.transportNonLinearSolver.LinearSolver = BackslashSolverAD();

reorder = pack(modelreorder, 'reorder', nls, 'Reordering');

%%

% Cartesian coarse grid
G = computeCellDimensions2(G);
G_cart = cartGrid([100, 100]);
p_cart = partitionUI(G_cart, [20, 20]);
p_cart = sampleFromBox(G, reshape(p_cart, G_cart.cartDims));
GC = generateCoarseGrid(G, p_cart);
GC = coarsenGeometry(GC);
GC = coarsenCellDimensions(GC);

mappings = getRefinementMappings(GC, GC, G, [1, GC.cells.num]);
GC        = generateCoarseGrid(GF, mappings.newPartition);
GC        = coarsenGeometry(GC);
GC.cells.refined = mappings.refined;

% [G, map1] = refineGrid(GC, GC, GF, [1, GC.cells.num]');
CG = coarsenCellDimensions(GC);

GC.cells.ghost = false(GC.cells.num,1);
G.cells.ghost = false(G.cells.num,1);

modeladapt = AdaptiveSequentialPressureTransportModel(modelseq.pressureModel, modelseq.transportModel, GC);

%%

problems = {fim, seq, reorder};

%% Simulate problems

runIx = 3;
for pNo = runIx
    [ok, status] = simulatePackedProblem(problems{pNo});
end

%%

setup = problems{3}.SimulatorSetup;
[ws, st, rep] = simulateScheduleAD(setup.state0, setup.model, setup.schedule);
