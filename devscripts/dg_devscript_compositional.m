mrstModule add msrsb compositional dg gasinjection ad-blackoil ...
    blackoil-sequential ad-core ad-props coarsegrid mrst-gui vista vemmech
mrstVerbose on

%%

caseName = 'simple_1d_wat';
% caseName = 'simple_3ph';
caseName = 'immiscible_denis'
[state0, modelFI, schedule, CG] = setupCompositionalPaperCase(caseName);
modelFI = NaturalVariablesCompositionalModel(modelFI.G, modelFI.rock, modelFI.fluid, modelFI.EOSModel.fluid, 'water', modelFI.water);

%%

G     = modelFI.G;
rock  = modelFI.rock;
fluid = modelFI.fluid;
compFluid = modelFI.EOSModel.fluid;
G = computeCellDimensions2(G);
G.cells.ghost = false(G.cells.num,1);
args = {G, rock, fluid, modelFI.EOSModel.fluid, 'water', modelFI.water, 'nonlineartolerance', 1e-4};
useNat = true;
if useNat
    model = NaturalVariablesCompositionalModel(args{:});
    pmodel = PressureNaturalVariablesModel(args{:});
    tmodel = TransportNaturalVariablesModel(args{:});
else
    model = OverallCompositionCompositionalModel(args{:});
    pmodel = PressureOverallCompositionModel(args{:});
    tmodel = TransportOverallCompositionModel(args{:});
end
modelSeq = SequentialPressureTransportModel(pmodel, tmodel);

%%

ix = 1:numel(schedule.step.val);
subschedule = schedule;
subschedule.step.val     = subschedule.step.val(ix);
subschedule.step.control = subschedule.step.control(ix);
[wsFV, stFV, rep] = simulateScheduleAD(state0, modelSeq, subschedule);

%%

close all
plotToolbar(modelSeq.transportModel.G, stFV);
plotWellSols(wsFV);

%%

degree = 1;
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
modelDG.transportModel = TransportNaturalVariablesModelDG(G, rock, fluid, compFluid, ...
                                   'disc'    , disc    , ...
                                   'dsMaxAbs', 0.2     , ...
                                   'water', modelFI.water, ...
                                   'nonlinearTolerance', 1e-4, ...
                                   'useIncTolComposition', false);
% modelDG.pressureModel.disc = disc;

state0 = assignDofFromState(disc, state0);
[wsDG, stDG, rep] = simulateScheduleAD(state0, modelDG, subschedule);

%%
figure
plotToolbar(modelDG.transportModel.G, stDG);
plotWellSols({wsFV, wsDG});