mrstModule add dg vem vemmech ad-props ad-core ad-blackoil ...
    blackoil-sequential reorder matlab_bgl upr coarsegrid mrst-gui ...
    mrst-utils weno

mrstVerbose on

%% Common stuff

baseName = 'qfs-reorder';
location = 'home';
% location = 'work';
switch location
    case 'work'
        dataDir  = fullfile('/media/strene/806AB4786AB46C92/mrst-dg/rsc-2019', baseName);
    case 'home'
        dataDir  = fullfile('/home/oysteskl/Documents/phd/repositories/mrst-dg/simulation-output/rsc-2019', baseName);
end

%% make PEBI grid

n  = 50;
l = 1000;
G = pebiGrid(l/n, [l,l]);
close all
plotGrid(G)
axis equal tight
G = computeGeometry(G);
G = computeVEMGeometry(G);
G = computeCellDimensions2(G);

%% Set up model

gravity reset on
gravity([0, -9.81]);

rock = makeRock(G, 100*milli*darcy, 1);

fluid = initSimpleADIFluid('phases', 'WO'                        , ...
                           'rho'   , [100, 1000]*kilogram/meter^3, ...
                           'mu'    , [0.5, 0.5]*centi*poise      , ...
                           'n'     , [2, 2]                      );

% FI model
modelFI = TwoPhaseOilWaterModel(G, rock, fluid, 'gravity', g);
% Seq model
modelSI = getSequentialModelFromFI(modelFI);
% Reordering model
modelReorder = getSequentialModelFromFI(modelFI);
modelReorder.transportModel = ReorderingModel(modelReorder.transportModel);
modelReorder.transportModel.chunkSize = 5;
modelReorder.transportModel.buffer = 0;

% Set up discretization
[mt, ot, jt] = deal(0);
disc = DGDiscretization(modelSI.transportModel, ...
               'degree'               , 0          , ...
               'basis'                , 'legendre' , ...
               'useUnstructCubature'  , true       , ...
               'jumpTolerance'        , jt         , ...
               'outTolerance'         , ot         , ...
               'outLimiter'           , 'kill'     , ...
               'meanTolerance'        , mt         );
modelDG = getSequentialModelFromFI(modelFI);
modelDG.transportModel = TransportOilWaterModelDG(G, rock, fluid, 'disc', disc, 'gravity', g);
modelDG.transportModel = ReorderingModelDG(modelDG.transportModel);
modelDG.transportModel.chunkSize = 50;

%% Set up schedule

time = 2*year;
rate = 1.5*sum(poreVolume(G, rock))/time;
xw = [0,0; l,l];

sW     = 0.0;

W = [];
d = pdist2(G.cells.centroids, xw);

[~, ix] = min(d(:,1));
W = addWell(W, G, rock, ix, 'type', 'rate', 'val', rate    , 'comp_i', [1,0]);
[~, ix] = min(d(:,2));
W = addWell(W, G, rock, ix, 'type', 'bhp' , 'val', 50*barsa, 'comp_i', [1,0]);

dt    = 7*day;
dtvec = rampupTimesteps(time, dt, 0);
schedule = simpleSchedule(dtvec, 'W', W);

state0 = initResSol(G, 100*barsa, [sW,1-sW]);
state0 = assignDofFromState(disc, state0);

%%

reorder = packSimulationProblem(state0, modelDG, schedule, baseName, ...
                                                 'Directory', dataDir, ...
                                                 'Name'     , 'reorder'  );
seq = packSimulationProblem(state0, modelSI, schedule, baseName, ...
                                                'Directory', dataDir, ...
                                                'Name'     , 'seq'  );
problems = {reorder, seq};

%%

runIx = 1:numel(problems);
% runIx = 2;
for pNo = runIx
    [ok, status] = simulatePackedProblem(problems{pNo});
end

%%
                
setup = problems{1}.SimulatorSetup;
[ws, st, rep] = simulateScheduleAD(setup.state0, setup.model, setup.schedule);
