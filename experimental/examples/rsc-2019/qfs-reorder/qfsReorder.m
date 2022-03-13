mrstModule add dg vem vemmech ad-props ad-core ad-blackoil ...
    blackoil-sequential reorder matlab_bgl upr coarsegrid mrst-gui ...
    mrst-utils weno

mrstVerbose on

%% Common stuff

baseName = 'qfs-reorder';
% location = 'home';
location = 'work';
switch location
    case 'work'
        dataDir  = fullfile('/media/strene/806AB4786AB46C92/mrst-dg/rsc-2019', baseName);
    case 'home'
        dataDir  = fullfile('/home/oysteskl/Documents/phd/repositories/mrst-dg/simulation-output/rsc-2019', baseName);
end

pack = @(state0, model, schedule, name) ...
    packSimulationProblem(state0, model, schedule, baseName, ...
                                       'Name'     , name   , ...
                                       'Directory', dataDir);
                                            

%% make PEBI grid

n  = 8;
l = 1000;
G = pebiGrid(l/n, [l,l]);
close all
plotGrid(G)
axis equal tight
G = computeGeometry(G);
G = computeVEMGeometry(G);
G = computeCellDimensions2(G);

%% make rock/fluid, set up schedule

rock = makeRock(G, 100*milli*darcy, 1);

fluid = initSimpleADIFluid('phases', 'WO'                        , ...
                           'rho'   , [1, 1000]*kilogram/meter^3, ...
                           'mu'    , [0.5, 0.5]*centi*poise      , ...
                           'n'     , [2, 2]                      );

time = 2*year;
rate = 1.5*sum(poreVolume(G, rock))/time;
xw = [0,0; l,l];

sW     = 0.0;

W = [];
d = pdist2(G.cells.centroids, xw);

[~, ix] = min(d(:,1));
W = addWell(W, G, rock, ix, 'type', 'rate', 'val', rate    , 'comp_i', [1,0]);
[~, ix] = min(d(:,2));
W = addWell(W, G, rock, ix, 'type', 'bhp' , 'val', 50*barsa, 'comp_i', [1,0]);

dt    = 7*day;
dtvec = rampupTimesteps(time, dt, 0);
schedule = simpleSchedule(dtvec, 'W', W);

state0 = initResSol(G, 100*barsa, [sW,1-sW]);

%% Set up models

chunk = 1;

% gravity reset off

% % FI model
% modelFI = TwoPhaseOilWaterModel(G, rock, fluid);
% % Seq model
% modelSI = getSequentialModelFromFI(modelFI);
% % Set up discretization
% [mt, ot, jt] = deal(0);
% disc = DGDiscretization(modelSI.transportModel, ...
%                'degree'               , 0          , ...
%                'basis'                , 'legendre' , ...
%                'useUnstructCubature'  , true       , ...
%                'jumpTolerance'        , jt         , ...
%                'outTolerance'         , ot         , ...
%                'outLimiter'           , 'kill'     , ...
%                'meanTolerance'        , mt         );
% state0 = assignDofFromState(disc, state0);
% modelDG = getSequentialModelFromFI(modelFI);
% modelDG.transportModel = TransportOilWaterModelDG(G, rock, fluid, 'disc', disc);
% % Reordering model
% modelDG.transportModel = ReorderingModelDG(modelDG.transportModel);
% modelDG.transportModel.chunkSize = chunk;
% modelDG.transportModel.plotProgress = false;
% modelDG.transportModel.nonlinearTolerance = 1e-3;
% 
% reorder = pack(state0, modelDG, schedule, 'reorder');
% 
gravity reset on
gravity([0, -100]);
% 
% % FI model
% modelFI = TwoPhaseOilWaterModel(G, rock, fluid);
% % Seq model
% modelSI = getSequentialModelFromFI(modelFI);
% % Set up discretization
% [mt, ot, jt] = deal(0);
% disc = DGDiscretization(modelSI.transportModel, ...
%                'degree'               , 0          , ...
%                'basis'                , 'legendre' , ...
%                'useUnstructCubature'  , true       , ...
%                'jumpTolerance'        , jt         , ...
%                'outTolerance'         , ot         , ...
%                'outLimiter'           , 'kill'     , ...
%                'meanTolerance'        , mt         );
modelDG = getSequentialModelFromFI(modelFI);
modelDG.transportModel = TransportOilWaterModel(G, rock, fluid, 'disc', disc);
% Reordering model
modelDG.transportModel = ReorderingModel(modelDG.transportModel);
modelDG.transportModel.chunkSize = chunk;
modelDG.pressureModel.extraStateOutput = true;

reorderGravity = pack(state0, modelDG, schedule, 'gravity');

problems = {reorder, reorderGravity};

%%

runIx = 1:numel(problems);
% runIx = 2;
for pNo = runIx
    [ok, status] = simulatePackedProblem(problems{pNo});
end

%%

setup = problems{2}.SimulatorSetup;
[ws, stFI, rep] = simulateScheduleAD(setup.state0, setup.model, setup.schedule);

%%

[ws, stFI, rep] = simulateScheduleAD(state0, modelSI, schedule);

%% Copyright Notice
%
% <html>
% <p><font size="-1">
% Copyright 2009-2022 SINTEF Digital, Mathematics & Cybernetics.
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
