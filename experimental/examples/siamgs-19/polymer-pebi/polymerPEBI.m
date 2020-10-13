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

[G, rock, state0, schedule] = setupPolymerPEBI(fluid, 'n', 20, 'dt', 20*day, 'nRampup', 8);

%% Fully implicit problem

modelFI   = OilWaterPolymerModel(G, rock, fluid);
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

xc = G.cells.centroids;
xw = xc(vertcat(schedule.control(1).W.cells));
d = pdist2(xc, xw);
c = any(d < 50,2);
modelDG = getSequentialModelFromFI(modelFI);
modelDG.transportModel ...
    = TransportOilWaterPolymerModelDG(G, rock, fluid, ...
                                'dsMaxAbs'     , 0.2, ...
                                'degree'       , 1  , ...
                                'meanTolerance', mt , ...
                                'outTolerance' , ot , ...
                                'jumpTolerance', Inf, ...
                                'useUnstructCubature', true);
modelDG.transportModel.disc.degree0 = ones(G.cells.num,1);
modelDG.transportModel.disc.degree0(c) = 0;
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

[GRef, rockRef, state0Ref, scheduleRef] = setupPolymerPEBI(fluid, 'n', 40, 'dt', 20*day, 'nRampup', 8);
modelFIRef = OilWaterPolymerModel(GRef, rockRef, fluid);
problemFIRef = pack('fi-ref', state0Ref, modelFIRef, scheduleRef);

%% Simulate refined problems

simulatePackedProblem(problemFIRef);

%% Simulate dG(0) problems

simulatePackedProblem(problemDG0);

%% Simulate DG problems

simulatePackedProblem(problemDG1);

%% Load dG(0) results

[wsDG0, stDG0, repDG0] = getPackedSimulatorOutput(problemDG0);

%% Load dG(1) results

[wsDG1, stDG1, repDG1] = getPackedSimulatorOutput(problemDG1);

%% Load WENO results

[wsWENO, stWENO, repWENO] = getPackedSimulatorOutput(problemWENO);

%% Load reference results

[wsFIRef, stFIRef, repFIRef] = getPackedSimulatorOutput(problemFIRef);

%%

close all

setup = problemSI.SimulatorSetup;
st    = stSI;

plotToolbar(setup.model.G, st, 'plot1d', true);

%%

close all
plotWellSols([wsSI, wsDG1], schedule.step.val)

%% Debug

setup = problemDG1.SimulatorSetup;
% setup = problemWENO.SimulatorSetup;
% setup = problemFV{1}.SimulatorSetup;
% setup = struct('model', model, 'state0', state0, 'schedule', schedule);
[ws, st, rep] = simulateScheduleAD(setup.state0, setup.model, setup.schedule);

%% Copyright Notice
%
% <html>
% <p><font size="-1">
% Copyright 2009-2020 SINTEF Digital, Mathematics & Cybernetics.
% </font></p>
% <p><font size="-1">
% This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).
% </font></p>
% <p><font size="-1">
% MRST is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% </font></p>
% <p><font size="-1">
% MRST is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% </font></p>
% <p><font size="-1">
% You should have received a copy of the GNU General Public License
% along with MRST.  If not, see
% <a href="http://www.gnu.org/licenses/">http://www.gnu.org/licenses</a>.
% </font></p>
% </html>
