mrstModule add dg vem vemmech ad-props ad-core ad-blackoil upr...
    blackoil-sequential mrst-gui reorder matlab_bgl ad-eor trust-region ...
    weno deckformat
mrstVerbose on
gravity reset off

%%

baseName = 'polymer-pebi';
dataDir  = fullfile(mrstPath('dg'), 'examples', 'siamgs-19', 'data', baseName);
pack = @(name, state0, model, schedule) ...
            packSimulationProblem(state0, model, schedule, baseName, ...
                                               'Directory', dataDir, ...
                                               'Name'     , name   );

%% Set up fluid

if 0
    fluid = initSimpleADIFluid('phases', 'WO'                      , ...
                               'rho'   , [1000, 1]*kilogram/meter^3, ...
                               'mu'    , [0.5 , 1]*centi*poise     , ...
                               'n'     , [1   , 1]                 );
    fluid = addSimplePolymerProperties(fluid, 'cmax', 5);
else
    fn   = fullfile(mrstPath('dg'), 'examples', 'siamgs-19', 'data', 'POLYMERDATA.DATA');
    deck = readEclipseDeck(fn);
    %The deck variables are converted to SI units.
    deck = convertDeckUnits(deck);
    % Setup the fluid structure
    fluid = initDeckADIFluid(deck);
    fluid.effads = @(c, cmax) effads(c, cmax, fluid);
end

%%

[G, rock, state0, schedule] = setupPolymerPEBI(fluid, 'n', 10, 'dt', 20*day, 'nRampup', 8);

%% Fully implicit problem

modelFI   = OilWaterPolymerModel(G, rock, fluid);
problemFI = pack('fi', state0, modelFI, schedule);

%% Run FI problem

simulatePackedProblem(problemFI);

%% Sequential problem

modelSI   = getSequentialModelFromFI(modelFI);
problemSI = pack('si', state0, modelSI, schedule);

%% dG(0) problem

[mt, ot] = deal(1e-5);
jt = 0.1;
modelDG = getSequentialModelFromFI(modelFI);
modelDG.transportModel ...
    = TransportOilWaterPolymerModelDG(G, rock, fluid, ...
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

modelDG = modelSI;
modelDG.transportModel ...
    = TransportOilWaterPolymerModelDG(G, rock, fluid, ...
                                'dsMaxAbs'     , 0.2, ...
                                'degree'       , 1  , ...
                                'meanTolerance', mt , ...
                                'outTolerance' , ot , ...
                                'jumpTolerance', jt , ...
                                'useUnstructCubature', true);
state0     = assignDofFromState(modelDG.transportModel.disc, state0);
problemDG1 = pack('dg-1', state0, modelDG, schedule);

%% Run dG(1) problem

simulatePackedProblem(problemDG1);

%% WENO problem

if 1
    modelWENO = HigherOrderOilWaterPolymerModel(G, rock, fluid);
%     modelWENO = HigherOrderOilWaterModel(G, rock, fluid);
    weno      = WENODiscretization(modelWENO, G.griddim, ...
                                 'includeBoundary'     , true, ...
                                 'interpolateReference', true);
    modelWENO.FluxDiscretization.saturationDiscretization = weno;
    modelWENO.FluxDiscretization.relPermDiscretization = weno;
    modelWENO.FluxDiscretization.discritizeRelPerm     = false;
    modelWENO.FluxDiscretization.discritizeViscosity   = false;
    modelWENO.dpMaxRel = 0.1;
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

[GRef, rockRef, state0Ref, scheduleRef] = setupPolymerPEBI(fluid, 'n', 30, 'dt', 20*day, 'nRampup', 8);
modelFIRef = OilWaterPolymerModel(GRef, rockRef, fluid);
problemFIRef = pack('fi-ref', state0Ref, modelFIRef, scheduleRef);

%% Simulate refined problems

simulatePackedProblem(problemFIRef);

%% Simulate dG(0) problems

simulatePackedProblem(problemDG0);

%% Simulate DG problems

simulatePackedProblem(problemDG1);

%% Load reference results

[wsRef, stRef, repRef] = getPackedSimulatorOutput(problemRef);

%% Load FI results

[wsFI, stFI, repFI] = getPackedSimulatorOutput(problemFI);

%% Load FV results

[wsSI, stSI, repSI] = getPackedSimulatorOutput(problemSI);

%% Load FV results

[wsDG0, stDG0, repDG0] = getPackedSimulatorOutput(problemDG0);

%% Load WENO results

[wsWENO, stWENO, repWENO] = getPackedSimulatorOutput(problemWENO);

%%

close all

setup = problemSI.SimulatorSetup;
st    = stSI;

plotToolbar(setup.model.G, st, 'plot1d', true);

%%

close all
plotWellSols([wsSI, wsDG(:,2)], schedule.step.val)

%% Debug

% setup = problemDG{1,2}.SimulatorSetup;
setup = problemWENO.SimulatorSetup;
% setup = problemFV{1}.SimulatorSetup;
% setup = struct('model', model, 'state0', state0, 'schedule', schedule);
[ws, st, rep] = simulateScheduleAD(setup.state0, setup.model, setup.schedule);