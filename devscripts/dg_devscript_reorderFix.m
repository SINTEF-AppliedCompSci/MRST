mrstModule add dg vem vemmech ad-props ad-core ad-blackoil ...
    blackoil-sequential vista mrst-gui reorder matlab_bgl weno
mrstVerbose on;

%%

gravity reset off

n = 10;
l = 1000;
G = computeGeometry(cartGrid([n,n], [l,l]*meter));
G = computeVEMGeometry(G);
G = computeCellDimensions(G);
[G.cells.equal, G.faces.equal] = deal(true);

rock = makeRock(G, 100*milli*darcy, 1);
fluid = initSimpleADIFluid('phases', 'WO'                   , ...
                           'rho'   , [1000, 1]*kilogram/meter^3, ...
                           'mu'    , [0.5, 0.5]*centi*poise     , ...
                           'n'     , [1, 1]                 );

modelfi = TwoPhaseOilWaterModel(G, rock, fluid);
modelFV = getSequentialModelFromFI(modelfi);
modelDG = modelFV;

%%

time = 2*year;
rate = 2*sum(poreVolume(G, rock))/time;
W = [];
W = addWell(W, G, rock, 1          , 'type', 'rate', 'val', rate    , 'comp_i', [1,0]);
W = addWell(W, G, rock, G.cells.num, 'type', 'bhp' , 'val', 50*barsa, 'comp_i', [1,0]);

if 0
    src = addSource([], 20, rate/3, 'sat', [1,0]);
else
    src = [];
end

dt    = 30*day;
dtvec = rampupTimesteps(time, dt, 0);
schedule = simpleSchedule(dtvec, 'W', W, 'src', src);

sW     = 0.0;
state0 = initResSol(G, 100*barsa, [sW,1-sW]);
[state0.posFlux, state0.negFlux] = deal(false);
% state0.cells = (1:G.cells.num);

%%

degree = 0:1;


% degree = [4];
% degree = [0,1,2];
[jt, ot, mt] = deal(Inf);
% 
% jt = Inf;
mt = 0;
ot = 1e-3;
% ot = 0.2;


% [wsDG, statesDG, wsDGreorder, statesDGreorder, disc] = deal(cell(numel(degree),1));
for dNo = 1%:numel(degree)
    disc{dNo} = DGDiscretization(modelDG.transportModel                   , ...
                                    'degree'               , degree(dNo)  , ...
                                    'basis'                , 'legendre'   , ...
                                    'useUnstructCubature'  , true        , ...
                                    'jumpTolerance'        , jt           , ...
                                    'outTolerance'         , ot           , ...
                                    'outLimiter'           , 'orderReduce', ...
                                    'meanTolerance'        , mt           , ...
                                    'limitAfterConvergence', false         , ...
                                    'plotLimiterProgress'  , false        );
    modelDG.transportModel = TransportOilWaterModelDG(G, rock, fluid, ...
                                       'disc'    , disc{dNo}        , ...
                                       'dsMaxAbs', 0.2, ...
                                       'nonlinearTolerance', 1e-3);
    modelDG.transportModel.conserveOil   = true;
    modelDG.transportModel.conserveWater = false;
    state0 = assignDofFromState(disc{dNo}, state0);
    [wsDG{dNo}, statesDG{dNo}, rep] = simulateScheduleAD(state0, modelDG, schedule); 
    
    modelDG.transportModel = ReorderingModelDG(modelDG.transportModel);
    modelDG.transportModel.chunkSize = 5;

    [wsDGreorder{dNo}, statesDGreorder{dNo}, rep] = simulateScheduleAD(state0, modelDG, schedule);
end

%%

m = TwoPhaseOilWaterModel(G, rock, fluid);
m = getSequentialModelFromFI(m);
[ws, st, r] = simulateScheduleAD(state0, m, schedule)

m.transportModel = ReorderingModel(m.transportModel);
m.transportModel.chunkSize = 5;
[wsr, str, rr] = simulateScheduleAD(state0, m, schedule);

%%

sdd = cellfun(@(s1,s2) compareStates(s1, s2) , st, str, 'unif', false);

%%

sd = cell(numel(degree),1);
for dNo = 1
    sd{dNo} = cellfun(@(s1,s2) compareStates(s1, s2), statesDG{dNo}, statesDGreorder{dNo}, 'unif', false);
end


%%

plotWellSols({wsDG{1}, wsDGreorder{1}});


