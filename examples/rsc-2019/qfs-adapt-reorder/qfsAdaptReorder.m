mrstModule add dg vem vemmech ad-props ad-core ad-blackoil ...
    blackoil-sequential gasinjection reorder matlab_bgl upr coarsegrid ...
    mrst-gui mrst-utils msrsb incomp weno

mrstVerbose on

%% Common stuff

baseName = 'qfs-dadapt-reorder';
% location = 'home';
location = 'work';
switch location
    case 'work'
        dataDir  = fullfile('/media/strene/806AB4786AB46C92/mrst-dg/rsc-2019', baseName);
    case 'home'
        dataDir  = fullfile('/home/oysteskl/Documents/phd/repositories/mrst-dg/simulation-output/rsc-2019', baseName);
end

pack = @(state0, model, schedule, name) ...
    packSimulationProblem(state0, model, schedule, baseName, ...
                                       'Name'     , name   , ...
                                       'Directory', dataDir);

%% Make PEBI grid

n  = 50;
l = 1000;
GF = pebiGrid(l/n, [l,l]);
close all
plotGrid(GF)
axis equal tight

GF = computeGeometry(GF);
GF = computeCellDimensions2(GF);

%% Make coarse grid

if 0
    rng(0);
    K = gaussianField([100, 100]);
    K = K - min(K(:));
    K = 100*milli*darcy*K./mean(K(:));
    K = max(K, 1*milli*darcy);
    K = sampleFromBox(GF, K);
else
    K = 100*milli*darcy;
end
rockF = makeRock(GF, K, 1);
T = computeTrans(GF, rockF);

p = partitionMETIS(GF, T, 300);
GC = generateCoarseGrid(GF, p);
GC = coarsenGeometry(GC);
GC = coarsenCellDimensions(GC);

wc = pdist2(GC.cells.centroids,[0,0;l,l]);
[~, wc] = min(wc);

mappings = getRefinementMappings(GC, GC, GF, wc);
G        = generateCoarseGrid(GF, mappings.newPartition);
G        = coarsenGeometry(G);
G.cells.refined = mappings.refined;
[G.cells.equal, G.faces.equal] = deal(false);
[GF.cells.equal, GF.faces.equal] = deal(false);

G = coarsenCellDimensions(G);

G.cells.ghost = false(G.cells.num,1);
GF.cells.ghost = false(GF.cells.num,1);

plotGrid(G)

%% Set up models

gravity reset off

rock  = makeRock(G, 100*milli*darcy, 1);
fluid = initSimpleADIFluid('phases', 'WO'                      , ...
                           'rho'   , [1000, 1]*kilogram/meter^3, ...
                           'mu'    , [0.5, 0.5]*centi*poise    , ...
                           'n'     , [1, 1]                    );

modelFIF = TwoPhaseOilWaterModel(GF, rockF, fluid);
modelSIF = getSequentialModelFromFI(modelFIF);

modelFI  = upscaleModelTPFA(modelFIF, G.partition);
modelSI  = getSequentialModelFromFI(modelFI);

[jt, ot, mt] = deal(Inf);
mt = 1e-8;
ot = 1e-8;
jt = 0.1;
% mt = 0; ot = 0;
degree = 0;
disc = DGDiscretization(modelSIF.transportModel               , ...
                        'degree'               , degree       , ...
                        'basis'                , 'legendre'   , ...
                        'useUnstructCubature'  , true         , ...
                        'jumpTolerance'        , jt           , ...
                        'outTolerance'         , ot           , ...
                        'outLimiter'           , 'orderReduce', ...
                        'meanTolerance'        , mt           , ...
                        'limitAfterConvergence', false        , ...
                        'plotLimiterProgress'  , false        );
transportModelDG = TransportOilWaterModelDG(GF, rockF, fluid, ...
                                   'disc'    , disc        , ...
                                   'dsMaxAbs', 0.2, ...
                                   'nonlinearTolerance', 1e-3);

tm = ReorderingModelDG(transportModelDG);
modelASIDGreorder = AdaptiveSequentialPressureTransportModel(modelSIF.pressureModel, tm, G);
modelASIDGreorder.transportModel.chunkSize = 50;
modelASIDGreorder.transportModel.parent.extraStateOutput = true;

%% Set up schedule

time  = 2*year;
dt    = 20*day;
dtvec = rampupTimesteps2(time, dt, 0);
rate  = 1.5*sum(poreVolume(GF, rockF))/time;

xw = [0,0; l,l];

% Coarse wells
W = [];
d = pdist2(G.cells.centroids, xw);
[~, ix] = min(d(:,1));
W = addWell(W, G, rock, ix, 'type', 'rate', 'val', rate    , 'comp_i', [1,0]);
[~, ix] = min(d(:,2));
W = addWell(W, G, rock, ix, 'type', 'bhp' , 'val', 50*barsa, 'comp_i', [1,0]);
% Coarse schedule
schedule = simpleSchedule(dtvec, 'W', W);
% Coarse state0
sW     = 0.0;
state0 = initResSol(G, 100*barsa, [sW,1-sW]);
state0.bfactor = [fluid.bW(state0.pressure), fluid.bO(state0.pressure)];

% Fine wells
WF = [];
d = pdist2(GF.cells.centroids, xw);
[~, ix] = min(d(:,1));
WF = addWell(WF, GF, rockF, ix, 'type', 'rate', 'val', rate    , 'comp_i', [1,0]);
[~, ix] = min(d(:,2));
WF = addWell(WF, GF, rockF, ix, 'type', 'bhp' , 'val', 50*barsa, 'comp_i', [1,0]);
% Fine schedule
scheduleF = simpleSchedule(dtvec, 'W', WF);
% Fine state0
state0F = initResSol(GF, 100*barsa, [sW,1-sW]);
state0F.bfactor = [fluid.bW(state0F.pressure), fluid.bO(state0F.pressure)];

%%

state0F.transportState        = state0;
state0F.transportState.degree = degree*ones(G.cells.num, 1);
state0F.transportState.G      = G;
state0F.transportState.pv     = modelASIDGreorder.transportModel.parent.operators.pv;
state0F.transportState        = assignDofFromState(modelASIDGreorder.transportModel.parent.disc, state0F.transportState);
% state0F.transportState        = disc.updateDofPos(state0F.transportState);
state0F.transportModel        = modelASIDGreorder.transportModel;

modelASIDGreorder.transportModel.parent.extraStateOutput = true;
modelASIDGreorder.pressureModel.extraStateOutput = true;
modelASIDGreorder.transportModel.plotProgress = false;
modelASIDGreorder.storeGrids = true;
modelASIDGreorder.plotProgress = true;
modelASIDGreorder.computeCoarsePressure = true;
modelASIDGreorder.transportModel.nonlinearSolver.errorOnFailure = true;
modelASIDGreorder.transportModel.nonlinearSolver.continueOnFailure = true;
modelASIDGreorder.coarsenTol = 5e-2;
modelASIDGreorder.refineTol = 2e-2;

reorderdg0 = pack(state0F, modelASIDGreorder, scheduleF, 'reorder-dg0');

%%

problems = {reorderdg0};

%%

for pNo = 1:numel(problems)
    [ok, status] = simulatePackedProblem(problems{pNo});
end

%%

setup = problems{1}.SimulatorSetup;
[ws, st, rep] = simulateScheduleAD(setup.state0, setup.model, setup.schedule);

