mrstModule add dg vem vemmech ad-props ad-core ad-blackoil ...
    blackoil-sequential gasinjection mrst-gui reorder matlab_bgl coarsegrid
mrstVerbose on;

%%

n  = 4;
GF = computeGeometry(cartGrid([n,n]));
p  = partitionUI(GF, [2, 2]);
GC = generateCoarseGrid(GF, p);
GC = coarsenGeometry(GC);

%%

close all

G = refineGrid(GC, GF, [55,56]');
plotToolbar(GF, G.partition);
colormap(hsv)
axis equal tight
% G = computeVEMGeometry(G)

%%

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
% state0.cells = (1:G.cells.num);

%%

degree = 0:5;


% degree = [4];
degree = [1,2];
[jt, ot, mt] = deal(Inf);
% 
% jt = Inf;
mt = 1e-6;
ot = 1e-3;
% ot = 0.2;


[wsDG, statesDG, disc] = deal(cell(numel(degree),1));
nls = NonLinearSolver('maxIterations', 25, 'useLinesearch', false);
for dNo = 1:numel(degree)
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
    
    modelDG.pressureModel = PressureOilWaterModelSemiDG(G, rock, fluid, ...
                                       'disc'    , disc{dNo}        );
    modelDG.transportNonLinearSolver = nls;

    state0 = assignDofFromState(modelDG.transportModel.disc, state0);
    [wsDG{dNo}, statesDG{dNo}, rep] = simulateScheduleAD(state0, modelDG, schedule);
end
