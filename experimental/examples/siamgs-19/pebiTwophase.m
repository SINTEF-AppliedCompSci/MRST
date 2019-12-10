mrstModule add dg vem vemmech ad-props ad-core ad-blackoil upr...
    blackoil-sequential mrst-gui reorder matlab_bgl ad-eor trust-region ...
    weno deckformat
mrstVerbose on
gravity reset on

%%

baseName = 'pebi-2ph';
dataDir  = fullfile(mrstPath('dg'), 'examples', 'siamgs-19', 'data', baseName);
pack = @(name, state0, model, schedule) ...
            packSimulationProblem(state0, model, schedule, baseName, ...
                                               'Directory', dataDir, ...
                                               'Name'     , name   );

%%

[G, rock, fluid, state0, schedule] = setupPEBItwoPhase('n', 20, 'dt', 20*day, 'nRampup', 8);

%% Fully implicit problem

modelFI   = TwoPhaseOilWaterModel(G, rock, fluid);
problemFI = pack('fi', state0, modelFI, schedule);

%% Run FI problem

simulatePackedProblem(problemFI);

%% Sequential problem

modelSI   = getSequentialModelFromFI(modelFI);
problemSI = pack('si', state0, modelSI, schedule);

%% Run SI problem

simulatePackedProblem(problemSI);

%% dG(0) problem

[mt, ot] = deal(1e-6);
% [mt, ot] = deal(0);
jt = 0.1;
modelDG = getSequentialModelFromFI(modelFI);
modelDG.transportModel ...
    = TransportBlackOilModelDG(G, rock, fluid, ...
                                'dsMaxAbs'     , 0.2, ...
                                'degree'       , 0  , ...
                                'meanTolerance', mt , ...
                                'outTolerance' , ot , ...
                                'useUnstructCubature', true);
state0     = assignDofFromState(modelDG.transportModel.disc, state0);
problemDG0 = pack('dg-0', state0, modelDG, schedule);

%% Run dG(0) problem

simulatePackedProblem(problemDG0);

%% dG(1) problem

modelDG = getSequentialModelFromFI(modelFI);
modelDG.transportModel ...
    = TransportBlackOilModelDG(G, rock, fluid, ...
                                'dsMaxAbs'     , 0.2, ...
                                'degree'       , 1  , ...
                                'meanTolerance', mt , ...
                                'outTolerance' , ot , ...
                                'jumpTolerance', jt, ...
                                'useUnstructCubature', true);
xc = G.cells.centroids;
xw = xc(vertcat(schedule.control(1).W.cells));
d = pdist2(xc, xw);
c = any(d < 30,2);
modelDG.transportModel.disc.degree0 = ones(G.cells.num,1);
modelDG.transportModel.disc.degree0(c) = 0;
state0     = assignDofFromState(modelDG.transportModel.disc, state0);
problemDG1 = pack('dg-1', state0, modelDG, schedule);

%% Run dG(1) problem

simulatePackedProblem(problemDG1);

%% WENO problem

if 1
    modelWENO = HigherOrderOilWaterModel(G, rock, fluid);
%     modelWENO = HigherOrderOilWaterModel(G, rock, fluid);
    weno      = WENODiscretization(modelWENO, G.griddim, ...
                                 'includeBoundary'     , true, ...
                                 'interpolateReference', true);
                             weno.implicit = true
                             weno.compactStencil = false
    modelWENO.FluxDiscretization.saturationDiscretization = weno;
    modelWENO.FluxDiscretization.relPermDiscretization    = weno;
%     modelWENO.FluxDiscretization.relPermDiscretization.implicit = true;
    modelWENO.FluxDiscretization.discritizeRelPerm     = false;
    modelWENO.FluxDiscretization.discritizeViscosity   = false;
    problemWENO = pack('weno', state0, modelWENO, schedule);
    
else
    foo = struct();
    [foo.cells, foo.faces, foo.sW, foo.c, foo.rate, foo.pressure] = deal([]);
    bc = struct();
    [bc.injection, bc.dirichlet] = deal(foo);
    c1 = G.faces.neighbors(G.faces.neighbors(:,1)==0,2);
    c2 = G.faces.neighbors(G.faces.neighbors(:,2)==0,1);
    G.cells.boundarycells = unique([c1; c2]);
    G = indexlist(G);
    modelWENO = PolymerHighOrderModelWell(G, rock, fluid, bc);
    modelWENO.do_high_order = true;
    modelWENO.reconstruction_type = 'weno';
    modelWENO.do_implicit_transport = true;
    modelWENO.do_primary = false;
    modelWENO.dpMaxRel = 0.1;
    
    scheduleWENO = schedule;
    [scheduleWENO.control(:).bc] = deal(bc);
    problemWENO = pack('weno', state0, modelWENO, scheduleWENO);
end

%% Simulate WENO

simulatePackedProblem(problemWENO);

%% Spatially refined problem

[GRef, rockRef, state0Ref, fluid, scheduleRef] = setupThreephasePEBI('n', 40, 'dt', 20*day, 'nRampup', 8);
modelFIRef = OilWaterPolymerModel(GRef, rockRef, fluid);
problemFIRef = pack('fi-ref', state0Ref, modelFIRef, scheduleRef);

%% Simulate refined problems

simulatePackedProblem(problemFIRef);

%%

close all

setup = problemDG0.SimulatorSetup;
st    = stSI;

plotToolbar(setup.model.G, st, 'plot1d', true);

%%

close all
plotWellSols([wsSI, wsDG1], schedule.step.val)

%% Debug

setup = problemDG0.SimulatorSetup;
% setup = problemWENO.SimulatorSetup;
% setup = problemFV{1}.SimulatorSetup;
% setup = struct('model', model, 'state0', state0, 'schedule', schedule);
[ws, st, rep] = simulateScheduleAD(setup.state0, setup.model, setup.schedule);
