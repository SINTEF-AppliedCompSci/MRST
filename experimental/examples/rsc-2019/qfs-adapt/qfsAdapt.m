mrstModule add dg vem vemmech ad-props ad-core ad-blackoil ...
    blackoil-sequential gasinjection reorder matlab_bgl upr coarsegrid ...
    mrst-gui mrst-utils weno

mrstVerbose on

%% Common stuff

dz = -3;
pw = @(G,W) plot3(G.cells.centroids([W.cells], 1), ...
                  G.cells.centroids([W.cells], 2), ...
                  G.cells.centroids([W.cells], 3) + dz, ...
                 'ok', 'markerSize', 8, 'markerFaceColor', 'w', 'lineWidth', 2);

baseName = 'qfs-adapt';
% location = 'home';
location = 'work';
switch location
    case 'work'
        dataDir  = fullfile('/media/strene/806AB4786AB46C92/mrst-dg/rsc-2019', baseName);
    case 'home'
        dataDir  = fullfile('/home/oysteskl/Documents/phd/repositories/mrst-dg/simulation-output/rsc-2019', baseName);
end

%% make PEBI grid

n  = 50;
l = 1000;
GF = pebiGrid(l/n, [l,l]);
close all
plotGrid(GF)
axis equal tight

%% Make coarse grid

GF = computeGeometry(GF);
GF = computeCellDimensions2(GF);

% Cartesian coarse grid
G_cart = cartGrid([100, 100]);
p_cart = partitionUI(G_cart, [10, 10]);
p_cart = sampleFromBox(GF, reshape(p_cart, G_cart.cartDims));
GC = generateCoarseGrid(GF, p_cart);
GC = coarsenGeometry(GC);
GC = coarsenCellDimensions(GC);

% Refine well cells
mappings = getRefinementMappings(GC, GC, GF, [1, GC.cells.num]);
G        = generateCoarseGrid(GF, mappings.newPartition);
G        = coarsenGeometry(G);
G.cells.refined = mappings.refined;
G = coarsenCellDimensions(G);
G.cells.ghost  = false(G.cells.num,1);
GF.cells.ghost = false(GF.cells.num,1);

%% Set up model

gravity reset off

rock  = makeRock(G, 100*milli*darcy, 1);
rockF = makeRock(GF, 100*milli*darcy, 1);

fluid = initSimpleADIFluid('phases', 'WO'                      , ...
                           'rho'   , [1000, 1]*kilogram/meter^3, ...
                           'mu'    , [0.5, 0.5]*centi*poise    , ...
                           'n'     , [1, 1]                    );

% FI model
modelFI = TwoPhaseOilWaterModel(GF, rockF, fluid);
% Seq model
modelSI = getSequentialModelFromFI(modelFI);
% ASI model
modelASI = AdaptiveSequentialPressureTransportModel(modelSI.pressureModel, modelSI.transportModel, G);
modelASI.transportModel.G = coarsenCellDimensions(modelASI.transportModel.G);
modelASI.refineTol  = 1*1e-2;
modelASI.coarsenTol = 3*1e-2;
% Coarse model
modelFIC = upscaleModelTPFA(modelFI, G);
modelSIC = getSequentialModelFromFI(modelFIC);

%% Set up schedule


time = 2*year;
dtvec = rampupTimesteps2(time, 7*day, 0);
rate = 1.5*sum(poreVolume(GF, rockF))/time;
xw = [0,0; l,l];

sW     = 0.0;
state0C = initResSol(G, 100*barsa, [sW,1-sW]);
state0C.bfactor = [fluid.bW(state0C.pressure), fluid.bO(state0C.pressure)];

WF = [];
d = pdist2(GF.cells.centroids, xw);

[~, ix] = min(d(:,1));
WF = addWell(WF, GF, rockF, ix, 'type', 'rate', 'val', rate    , 'comp_i', [1,0]);
[~, ix] = min(d(:,2));
WF = addWell(WF, GF, rockF, ix, 'type', 'bhp' , 'val', 50*barsa, 'comp_i', [1,0]);

schedule = simpleSchedule(dtvec, 'W', WF);

state0 = initResSol(GF, 100*barsa, [sW,1-sW]);
state0.bfactor = [fluid.bW(state0.pressure), fluid.bO(state0.pressure)];

WC = [];
d = pdist2(G.cells.centroids, xw);

[~, ix] = min(d(:,1));
WC = addWell(WC, G, rock, ix, 'type', 'rate', 'val', rate    , 'comp_i', [1,0]);
[~, ix] = min(d(:,2));
WC = addWell(WC, G, rock, ix, 'type', 'bhp' , 'val', 50*barsa, 'comp_i', [1,0]);
scheduleC = simpleSchedule(dtvec, 'W', WC);

%%

modelASI.storeGrids                     = true;
modelASI.plotProgress                   = true;
modelASI.pressureModel.extraStateOutput = true;

state0.transportState        = state0C;
% state0.transportState.degree = degree*ones(G.cells.num, 1);
state0.transportState.G      = G;
state0.transportState.pv     = modelASI.transportModel.operators.pv;
% state0.transportState        = assignDofFromState(discDG, state0.transportState);
% state0.transportState        = discDG.updateDofPos(state0.transportState);
state0.transportModel        = modelASI.transportModel;

%%

adaptive = packSimulationProblem(state0, modelASI, schedule, baseName, ...
                                                 'Directory', dataDir, ...
                                                 'Name'     , 'adapt'  );
seq = packSimulationProblem(state0, modelSI, schedule, baseName, ...
                                                'Directory', dataDir, ...
                                                'Name'     , 'seq'  );
coarse = packSimulationProblem(state0C, modelSIC, scheduleC, baseName, ...
                                                'Directory', dataDir, ...
                                                'Name'     , 'coarse'  );

problems = {adaptive, seq, coarse};

%%

runIx = 1:numel(problems);
for pNo = runIx
    [ok, status] = simulatePackedProblem(problems{pNo});
end

%%
                
setup = problems{3}.SimulatorSetup;
[ws, st, rep] = simulateScheduleAD(setup.state0, setup.model, setup.schedule);

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
