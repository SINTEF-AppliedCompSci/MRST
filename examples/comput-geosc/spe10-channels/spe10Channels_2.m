mrstModule add dg vem vemmech ad-props ad-core ad-blackoil blackoil-sequential gasinjection reorder matlab_bgl upr mrst-gui spe10
mrstVerbose on

%%

[state0, modelFI, schedule] = setupSPE10_AD('layers', 50);
G = modelFI.G;
G = computeVEMGeometry(G);
G = computeCellDimensions(G);
modelFI.G = G;
rock = modelFI.rock;
fluid = modelFI.fluid;
close all

time = 4*year;
rate = 0.2*sum(poreVolume(G, rock))/time;
bhp = 275*barsa;

W = [];
W = verticalWell(W, G, rock, 31, 1, [], ...
                 'type', 'rate', ...
                 'val', rate, ...
                 'comp_i', [1,0]);
W = verticalWell(W, G, rock, 35, 220, [], ...
                 'type', 'bhp', ...
                 'val', bhp, ...
                 'comp_i', [1,0]);
             
dt = 20*day;
dtvec = rampupTimesteps(time, dt);
schedule = simpleSchedule(dtvec, 'W', W);

%%

plotToolbar(modelFI.G, modelFI.rock);
plotWell(G, W);
axis equal tight

%%

modelFV = getSequentialModelFromFI(modelFI);
modelDG = modelFV;
[modelDG.transportModel.extraStateOutput, ...
 modelDG.pressureModel.extraStateOutput  ] = deal(true);

%%

ot = 1e-3;
jt = Inf;
mt = 0.0;

degree = [0,1,2,3];
dsMaxAbs = 0.1;

baseName = 'spe10-channels';
dataDir  = fullfile('/media/strene/806AB4786AB46C92/mrst-dg/comput-geosc', baseName);

runIx = 2:4;
[problems, reorderProblems] = deal(cell(numel(runIx), 1));
pNo = 1;
for dNo = runIx
    
    % Set up discretization
    disc   = DGDiscretization(modelDG.transportModel           , ...
                             'degree'             , degree(dNo), ...
                             'basis'              , 'legendre' , ...
                             'useUnstructCubature', true       ,  ...
                             'jumpTolerance'      , jt         , ...
                             'outTolerance'       , ot         , ...
                             'outLimiter'         , 'kill'     , ...
                             'meanTolerance'      , mt         , ...
                             'limitAfterConvergence', false    );
    % Set up transport model
    modelDG.transportModel = TransportOilWaterModelDG(G, rock, fluid, ...
                                'disc'              , disc    , ...
                                'nonlinearTolerance', 1e-3    , ...
                                'dsMaxAbs'          , dsMaxAbs);
    % Set up pressure model
%     modelDG.pressureModel = PressureOilWaterModelSemiDG(G, rock, fluid, ...
%                                 'disc'              , disc     , ...
%                                 'extraStateOutput'  , true);
    % Set up initial state
    state0 = assignDofFromState(modelDG.transportModel.disc, state0);
    % Add to problems
    problems{pNo} = packSimulationProblem(state0, modelDG, schedule, baseName, ...
                                          'Directory', dataDir, ...
                                          'Name'     , ['dg-', num2str(degree(dNo))']);
    % Make reordering model
    modelDGReorder = modelDG;
    modelDGReorder.transportModel ...
        = ReorderingModelDG_ghost(modelDGReorder.transportModel, ...
                                    'plotProgress'      , false, ...
                                    'plotAfterTimestep' , false, ...
                                    'plotAfterCellSolve', false, ...
                                    'chunkSize'         , 1    , ...
                                    'nonlinearTolerance', 1e-3 );
    modelDGReorder.transportModel.parent.nonlinearTolerance = 1e-3;
    modelDGReorder.transportModel.parent.extraStateOutput = true;
    % Add to problems
    reorderProblems{pNo} = packSimulationProblem(state0, modelDGReorder, schedule, baseName, ...
                                          'Directory', dataDir, ...
                                          'Name'     , ['reorder-dg-', num2str(degree(dNo))']);                 
    pNo = pNo+1;
                                      
end

%%

[okDG, statusDG] = simulatePackedProblem(problems);

%%

[okDGReorder, statusDGReorder] = simulatePackedProblem(reorderProblems);

%%