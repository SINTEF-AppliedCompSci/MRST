mrstModule add spe10 vem vemmech dg ad-core ad-props ad-blackoil ...
    blackoil-sequential vista incomp msrsb matlab_bgl coarsegrid ...
    mrst-gui reorder weno opm-data compositional linearsolvers

mrstVerbose on

%% Common stuff

dz = -3;
pw = @(G,W) plot3(G.cells.centroids([W.cells], 1), ...
                  G.cells.centroids([W.cells], 2), ...
                  G.cells.centroids([W.cells], 3) + dz, ...
                 'ok', 'markerSize', 8, 'markerFaceColor', 'w', 'lineWidth', 2);

baseName = 'norne';
% location = 'home';
location = 'work';
switch location
    case 'work'
        dataDir  = fullfile('/media/strene/806AB4786AB46C92/mrst-dg/rsc-2019', baseName);
    case 'home'
        dataDir  = fullfile('/home/oysteskl/Documents/phd/repositories/mrst-dg/simulation-output/rsc-2019', baseName);
end

pack = @(state0, model, schedule, name) ...
    packSimulationProblem(state0, model, schedule, baseName, 'Name', name, 'Directory', dataDir);

%% Get model

[state0, model, schedule] = getReorderingCase('norne_simple');

%% Make coarse grid

G     = model.G;
rock  = model.rock;
fluid = model.fluid;

G = computeCellDimensions2(G);

T = computeTrans(G, rock);
p = partitionMETIS(G, T, ceil(G.cells.num/10));

%% Make coarse grid

W  = schedule.control(1).W;
wc = vertcat(W.cells);
xw = G.cells.centroids(wc, :);
d  = pdist2(G.cells.centroids(:, 1:2), xw(:, 1:2));
refine = false(G.cells.num, 1);
refine = refine | any(d < 300*meter,2);
refine(wc) = true;
p(refine) = max(p) + find(refine);
p = compressPartition(p);

GC = generateCoarseGrid(G, p);
GC = coarsenGeometry(GC);
GC = coarsenCellDimensions(GC, 'useQhull', true);
GC = storeInteractionRegionCart(GC);

d  = pdist2(GC.cells.centroids(:, 1:2), xw(:, 1:2));
GC.cells.refined = any(d < 300*meter,2);


%% Sequential model

modelSeq = getSequential_local(model);
seq      = pack(state0, modelSeq, schedule, 'sequential');

%% 

modelSeq = getSequential_local(model);


modelAdaptive ...
    = AdaptiveSequentialPressureTransportModel(modelSeq.pressureModel, modelSeq.transportModel, GC);
modelAdaptive.refineTol = 1e-3;
modelAdaptive.coarsenTol = 0.5*1e-3;
tmodel = modelAdaptive.transportModel;


state0a = state0;
state0a = ensureDensitiesPresentCompositional(model, state0a);

coarsemodel = upscaleModelTPFA(model, GC);
state0C = upscaleState(tmodel, model, state0a);
% state0C = ensureDensitiesPresentCompositional(tmodel, state0C);

state0a.transportState = state0C;
state0a.transportState.pv = tmodel.operators.pv;
state0a.transportState.G = GC;
state0a.transportModel = tmodel;
state0a.wellSol = [];


modelAdaptive.computeCoarsePressure = false;
modelAdaptive.plotProgress = true;

adapt = pack(state0a, modelAdaptive, schedule, 'adaptive');

%%

problems = {seq, adapt};

%%

for pNo = 2
    [ok, status] = simulatePackedProblem(problems{pNo});
end

%%

setup = problems{2}.SimulatorSetup;
[ws, st, rep] = simulateScheduleAD(setup.state0, setup.model, setup.schedule);
