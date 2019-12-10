useWater = true;
caseName = 'verysimple';
p = 150*barsa;
bhp = 100*barsa;
pvi = 0.2;
nstep = 30;
time = 1*year;
n = 10;

[G, rock, state0, schedule, fluid, eos] = ...
    getTest1D_compositional(n, 'water', useWater, ...
                               'fluid', caseName, ...
                               'nkr', 2, ...
                               'nstep', nstep, ...
                               'pvi', pvi, ...
                               'pressure', p, ...
                               'bhp', bhp, ...
                               'perm', 0.5*darcy, ...
                               'time', time);
G = computeCellDimensions2(G);
G.cells.ghost = false(G.cells.num,1);

args = {G, rock, fluid, eos.fluid, 'water', useWater, 'useIncTolComposition', true};

pmodel = PressureNaturalVariablesModel(args{:});
tmodel = TransportNaturalVariablesModel(args{:});
modelSeq = SequentialPressureTransportModel(pmodel, tmodel);

ix = 1:30;
subschedule = schedule;
subschedule.step.val     = subschedule.step.val(ix);
subschedule.step.control = subschedule.step.control(ix);

%%

[wsFV, stFV, rep] = simulateScheduleAD(state0, modelSeq, subschedule);

%%

close all
plotToolbar(G, stFV, 'plot1d', true)

%%

degree = 0;
[jt, ot, mt] = deal(Inf);
mt = 0.0;
% pmodel = PressureNaturalVariablesModelSemiDG(args{:});
pmodel = PressureNaturalVariablesModel(args{:});
tmodel = TransportNaturalVariablesModelDG(args{:});
modelDG = SequentialPressureTransportModel(pmodel, tmodel);
disc = DGDiscretization(tmodel                                        , ...
                                'degree'               , degree       , ...
                                'basis'                , 'legendre'   , ...
                                'useUnstructCubature'  , true         , ...
                                'jumpTolerance'        , jt           , ...
                                'outTolerance'         , ot           , ...
                                'outLimiter'           , 'orderReduce', ...
                                'meanTolerance'        , mt           , ...
                                'limitAfterConvergence', false        , ...
                                'plotLimiterProgress'  , false        );
modelDG.transportModel = TransportNaturalVariablesModelDG(G, rock, fluid, eos.fluid, ...
                                   'disc'    , disc    , ...
                                   'dsMaxAbs', 0.1     , ...
                                   'nonlinearTolerance', 1e-3, ...
                                   'useIncTolComposition', true);
% modelDG.pressureModel.disc = disc;

state0 = assignDofFromState(disc, state0);
[wsDG, stDG, rep] = simulateScheduleAD(state0, modelDG, subschedule);

%%

close all
figure; plotToolbar(G, stFV, 'plot1d', true)
figure; plotToolbar(G, stDG, 'plot1d', true)

%%
plotWellSols({wsFV, wsDG})
