mrstModule add spe10 vem vemmech dg ad-core ad-props ad-blackoil ...
    blackoil-sequential gasinjection incomp msrsb matlab_bgl coarsegrid ...
    mrst-gui reorder weno

mrstVerbose on

%% Common stuff

dz = -3;
pw = @(G,W) plot3(G.cells.centroids([W.cells], 1), ...
                  G.cells.centroids([W.cells], 2), ...
                  G.cells.centroids([W.cells], 3) + dz, ...
                 'ok', 'markerSize', 8, 'markerFaceColor', 'w', 'lineWidth', 2);

baseName = 'spe10';
% location = 'home';
location = 'work';
switch location
    case 'work'
        dataDir  = fullfile('/media/strene/806AB4786AB46C92/mrst-dg/rsc-2019', baseName);
    case 'home'
        dataDir  = fullfile('/home/oysteskl/Documents/phd/repositories/mrst-dg/simulation-output/rsc-2019', baseName);
end
            
%% Set up base model

[state0, model, schedule]  = setupSPE10_AD('layers', 70);
G        = model.G;
rock     = model.rock;
fluid     = model.fluid;
% schedule = schedule;
schedule.step.val = rampupTimesteps(sum(schedule.step.val), 20*day, 8);
schedule.step.control = ones(numel(schedule.step.val), 1);
G        = computeCellDimensions2(G);
G.equal  = true;

%% Construct coarse grid

T = computeTrans(G, rock);
p = partitionMETIS(G, T, 800);
GC = generateCoarseGrid(G, p);
GC = coarsenGeometry(GC);
GC = coarsenCellDimensions(GC, 'useQhull', true);
GC = storeInteractionRegionCart(GC);

%% Make new wells

time = 2000*day;
bhp  = 50*barsa;
rate = 0.25*sum(poreVolume(G, rock))/time;

WF = [];
WF = verticalWell(WF, G, rock, 1, 1, []  , ...
                                   'type'  , 'rate', ...
                                   'val'   , rate  , ...
                                   'comp_i', [1,0] , ...
                                   'name'  , 'Inj' );
WF = verticalWell(WF, G, rock, 42, 220, [], ...
                                   'type'  , 'bhp' , ...
                                   'val'   , bhp   , ...
                                   'comp_i', [1,0] , ...
                                   'name'  , 'Prod');
schedule.control(1).W = WF;

%% Make coarse grid with refinement around wells

xw = G.cells.centroids([WF.cells], :);

d      = pdist2(GC.cells.centroids, xw);
refine = any(d < 50*meter,2);

mappings = getRefinementMappings(GC, GC, G, refine);
GC        = generateCoarseGrid(G, mappings.newPartition);
GC.cells.refined = mappings.refined;
GC        = coarsenGeometry(GC);
GC        = coarsenCellDimensions(GC, 'useQhull', true);
GC.equal  = false;

% %% Make well struct for the coarse grid
% 
% wc = GC.partition([WF.cells]);
W = WF;
for wNo = 1:numel(W)
    W(wNo).cells = GC.partition(W(wNo).cells);
end
% schedule = schedule;
% schedule.control(1).W  = W;

%% Plot setup

close all
plotToolbar(G, rock);
plotGrid(GC, 'facec', 'none', 'edgec', [1,1,1]*0.25, 'linew', 1.5);
hold on
pw(GC, W);
colormap(pink)
axis equal tight; ax = gca;

%% Set up models

degree = [0,1];
jt = Inf; ot = 1e-3; mt = 1e-3;

% Fully implicit fine model
G.cells.equal = true;
G.faces.equal = false;
GC.cells.equal = false;
GC.faces.equal = false;
modelFI = TwoPhaseOilWaterModel(G, rock, model.fluid);

[modelDG, modelDGreorder, modelDGadapt, modelDGadaptReorder] = deal(cell(numel(degree),1));

for dNo = 1:numel(degree)
    
    % Assign sequential models
    [modelDG{dNo}, modelDGreorder{dNo}, modelDGadapt{dNo}, modelDGadaptReorder{dNo}] ...
        = deal(getSequentialModelFromFI(modelFI));

    % Set up discretization
    disc = DGDiscretization(modelDG{dNo}.transportModel, ...
                   'degree'               , degree(dNo), ...
                   'basis'                , 'legendre' , ...
                   'useUnstructCubature'  , true       , ...
                   'jumpTolerance'        , jt         , ...
                   'outTolerance'         , ot         , ...
                   'outLimiter'           , 'kill'     , ...
                   'meanTolerance'        , mt         );
    % Transport model
    tmodel = TransportOilWaterModelDG(G, rock, model.fluid, ...
                                'disc'              , disc, ...
                                'dsMaxAbs'          , 0.1 , ...
                                'nonlinearTolerance', 1e-3);
                            
    % Reordering model
    tmodelReorder = ReorderingModelDG(tmodel, 'chunkSize', 1);
    
    % Regular model
    modelDG{dNo}.transportModel = tmodel;
    
    % Reordering model
    modelDGreorder{dNo}.transportModel = tmodelReorder;
    
    % Adaptive model
    modelDGadapt{dNo} ...
        = AdaptiveSequentialPressureTransportModel(modelDG{dNo}.pressureModel, tmodel, GC);
                    
    % Adaptive reordering mdoel
    modelDGadaptReorder{dNo} ...
        = AdaptiveSequentialPressureTransportModel(modelDG{dNo}.pressureModel, tmodelReorder, GC);
    
end

modelRef = getSequentialModelFromFI(modelFI);

%%

state0.bfactor = [fluid.bW(state0.pressure), fluid.bO(state0.pressure)];
state0C = initResSol(GC, state0.pressure(1), state0.s(1,:));
state0C.bfactor = [fluid.bW(state0C.pressure), fluid.bO(state0C.pressure)];

%% Set up problems

names = {'base', 'reorder', 'adapt', 'adapt-reorder'};
problems = cell(numel(degree), numel(modelDG{1}));
for dNo = 1:numel(degree)
    
    tmodel = modelDG{dNo}.transportModel;
    state0         = assignDofFromState(tmodel.disc, state0);
    state0.bfactor = [fluid.bW(state0.pressure), fluid.bO(state0.pressure)];
    
    tmodel = modelDGadapt{dNo}.transportModel;
    state0C = assignDofFromState(tmodel.disc, state0C);
    pv = tmodel.operators.pv;
    d = fullfile(dataDir, ['dg', num2str(degree(dNo))]);
    mdls = {modelDG{dNo}, modelDGreorder{dNo}, modelDGadapt{dNo}, modelDGadaptReorder{dNo}};
    
    for mNo = 1:numel(mdls)
        
        if isa(mdls{mNo}, 'AdaptiveSequentialPressureTransportModel')
            tmodel  = mdls{mNo}.transportModel;
            state0.transportState = state0C;
            state0.transportState.pv = pv;
            state0.transportState.G = GC;
            state0.transportModel = tmodel;
            if isa(mdls{mNo}.transportModel, 'ReorderingModel')
                mdls{mNo}.computeCoarsePressure = true;
            end
            mdls{mNo}.plotProgress = true;
        end
        mdls{dNo}.pressureModel.extraStateOutput = true;
        problems{dNo,mNo} = packSimulationProblem(state0, mdls{mNo}, schedule, baseName, 'Directory', d, 'name', names{mNo});
    end
end

%%

modelWENO = HigherOrderOilWaterModel(G, rock, fluid);
disc      = WENODiscretization(modelWENO, 3, 'includeBoundary', true, ...
                                             'interpolateReference', true);
modelWENO.FluxDiscretization.saturationDiscretization = disc;
modelWENO.FluxDiscretization.relPermDiscretization    = disc;
modelWENO.FluxDiscretization.discritizeRelPerm        = false;
modelWENO.FluxDiscretization.discritizeViscosity      = false;
modelWENO.FluxDiscretization.relPermDiscretization.implicit = true;
weno      = packSimulationProblem(state0, modelWENO, schedule, baseName, ...
                'Directory', dataDir, 'name', 'weno');

%%

mdlIx = 1;
for dNo = 2
    for mNo = mdlIx
        [ok, status] = simulatePackedProblem(problems{dNo,mNo});
    end
end

%%

[ok, status] = simulatePackedProblem(weno);

%% Debug

% setup = problems{2,3}.SimulatorSetup;
% setup.model.plotProgress = true;
setup = weno.SimulatorSetup;
[ws, st, rep] = simulateScheduleAD(setup.state0, setup.model, setup.schedule);


