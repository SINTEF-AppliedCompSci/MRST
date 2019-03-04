mrstModule add dg vem vemmech ad-props ad-core ad-blackoil blackoil-sequential gasinjection

%%

n = 10;
l = 1000;
G = computeGeometry(cartGrid([n,n,n], [l,l,l]*meter));
% G = computeGeometry(cartGrid([n,1,1], [l,10,10]*meter));
G = computeVEMGeometry(G);
G = computeCellDimensions2(G);
[G.cells.equal, G.faces.equal] = deal(true);

rock = makeRock(G, 100*milli*darcy, 1);
fluid = initSimpleADIFluid('phases', 'WO'                   , ...
                           'rho'   , [1, 1]*kilogram/meter^3, ...
                           'mu'    , [1, 1]*centi*poise     , ...
                           'n'     , [1, 1]                 );
% fluid.krW = @(s) fluid.krW(s).*(s>=0 & s<=1) + fluid.krW(1).*(s>1);
% fluid.krO = @(s) fluid.krO(s).*(s>=0 & s<=1) + fluid.krO(1).*(s>1);

modelfi = TwoPhaseOilWaterModel(G, rock, fluid);
modelFV = getSequentialModelFromFI(modelfi);
modelDG = modelFV;

% G = makeFaceCoords(G);

%%

time = 2*year;
rate = 1*sum(poreVolume(G, rock))/time;
W = [];
W = addWell(W, G, rock, 1          , 'type', 'rate', 'val', rate    , 'comp_i', [1,0]);
W = addWell(W, G, rock, G.cells.num, 'type', 'bhp' , 'val', 50*barsa, 'comp_i', [1,0]);

dt    = 30*day;
dtvec = rampupTimesteps(time, dt, 0);
schedule = simpleSchedule(dtvec, 'W', W);

sW     = 0.0;
state0 = initResSol(G, 100*barsa, [sW,1-sW]);
state0.cells = (1:G.cells.num)';

%%

degree = [1,2];
% degree = 2;
[wsDG, statesDG] = deal(cell(numel(degree),1));
for dNo = 1:numel(degree)
    disc    = DGDiscretization(modelDG.transportModel, ...
                                'degree'             , degree(dNo), ...
                                'basis'              , 'legendre' , ...
                                'useUnstructCubature', true      , ...
                                'jumpTolerance'      , Inf        , ...
                                'outTolerance'       , 0.0        , ...
                                'outLimiter'         , 'kill', ...
                                'meanTolerance'      , 0.0        );
%     disc.limitAfterConvergence = true;
    disc.limitAfterNewtonStep  = true;
    modelDG.transportModel = TransportOilWaterModelDG(G, rock, fluid  , ...
                                                      'disc'    , disc, ...x
                                                      'dsMaxAbs', 0.1 );    
    modelDG.pressureModel = PressureOilWaterModelSemiDG(G, rock, fluid  , ...
                                                      'disc'    , disc);    
    
    
%     modelDG.transportModel.AutoDiffBackend = DiagonalAutoDiffBackend();
    nls = NonLinearSolver('maxIterations', 10);
    modelDG.transportNonLinearSolver = nls;
    state0 = assignDofFromState(modelDG.transportModel.disc, state0);
    [wsDG{dNo}, statesDG{dNo}, rep] = simulateScheduleAD(state0, modelDG, schedule);
end

%%

[wsFV, statesFV, rep] = simulateScheduleAD(state0, modelFV, schedule);

%%

close all

nclr = 6;
for dNo = 1:numel(degree)
    figure
    plotToolbar(G, statesDG{dNo});
    
%     colormap(summer(nclr))
end

figure
plotToolbar(G, statesFV);
% colormap(summer(nclr))

%%
    

plotWellSols({wsDG{:}, wsFV})