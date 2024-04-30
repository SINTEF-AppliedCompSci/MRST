mrstModule add spe10 vem vemmech dg ad-core ad-props ad-blackoil ...
    blackoil-sequential incomp msrsb matlab_bgl coarsegrid mrst-gui ...
    reorder weno vista compositional

mrstVerbose on

gravity reset off

%% Common stuff

name = 'qfs_co2_small';
baseName = 'qfs-CO2';

location = 'work';
% location = 'home';
switch location
    case 'work'
        dataDir = fullfile('/media/strene/806AB4786AB46C92/mrst-dg/rsc-2019', baseName);
    case 'home'
        dataDir = fullfile(mrstPath('dg'), 'examples', 'rsc-2019', 'qfs-co2', 'sim-output');
end
[state0, model, schedule] = getReorderingCase(name);

pack = @(state0, model, name, desc) ...
    packSimulationProblem(state0, model, schedule, baseName, ...
                          'Name'           , name          , ...
                          'Directory'      , dataDir       , ...
                          'Description'    , desc          );

%% Fully-implicit

model.AutoDiffBackend = DiagonalAutoDiffBackend();
fluid = model.fluid;

fim = pack(state0, model, 'fim', 'Fully-implicit');

%% Sequential-implicit

modelseq = getSequential_local(model);
modelseq.transportNonLinearSolver.LinearSolver = BackslashSolverAD();
seq = pack(state0, modelseq, 'seq', 'Sequential');

%% Sequential-reordering

modelreorder = getSequential_local(model);
modelreorder.transportNonLinearSolver.LinearSolver = BackslashSolverAD();
modelreorder.pressureModel.extraStateOutput = true;
modelreorder.transportModel = ReorderingModel(modelreorder.transportModel);
modelreorder.transportModel.parent.extraStateOutput = true;

modelreorder.transportModel.chunkSize = 50;
modelreorder.transportModel.buffer = 0;

reorder = pack(state0, modelreorder, 'reorder', 'Reordering');

%% Adaptive

nls = NonLinearSolver();

modeladapt = getSequential_local(model);
modeladapt.pressureModel.extraStateOutput = true;

p = partitionCartGrid(model.G.cartDims, [15,15,1]);
G = model.G;
G = computeCellDimensions2(G);
GC = generateCoarseGrid(G, p);
GC = coarsenGeometry(GC);
GC = coarsenCellDimensions(GC);
m = getRefinementMappings(GC, GC, model.G, [1, GC.cells.num]);

GC = generateCoarseGrid(model.G, m.newPartition);
GC = coarsenGeometry(GC);
GC.cells.refined = m.refined;

coarsemodel = upscaleModelTPFA(modeladapt.transportModel, m.newPartition);
state0C = upscaleState(coarsemodel, modeladapt.transportModel, state0);
modeladapt = AdaptiveSequentialPressureTransportModel(modeladapt.pressureModel, modeladapt.transportModel, GC);
state0 = ensureDensitiesPresentCompositional(model, state0);
state0C.rho = repmat(state0.rho(1,:), GC.cells.num,1);
state0C.bfactor = state0C.rho./[fluid.rhoOS, fluid.rhoGS];
state0C.x = repmat(state0.x(1,:), GC.cells.num, 1);
state0C.y = repmat(state0.y(1,:), GC.cells.num, 1);
state0C.L = repmat(state0.L(1,:), GC.cells.num, 1);
state0C.K = repmat(state0.K(1,:), GC.cells.num, 1);
state0C.Z_V = repmat(state0.Z_V(1,:), GC.cells.num, 1);
state0C.Z_L = repmat(state0.Z_L(1,:), GC.cells.num, 1);

state0C.components = repmat(state0C.components(1,:), GC.cells.num, 1);
state0.transportState = state0C;
state0.transportState.G = GC;
state0.transportState.pv = coarsemodel.operators.pv;
state0.transportModel = modeladapt.transportModel;
modeladapt.plotProgress = true;
modeladapt.refineTol = 1e-2;
modeladapt.coarsenTol = 1e-2;

adapt = pack(state0, modeladapt, 'adapt', 'Adaptive');

%% Coarse 

coarsemodel = upscaleModelTPFA(model, m.newPartition);

state0C.bfactor = state0C.rho./[fluid.rhoOS, fluid.rhoGS];
state0C.x = repmat(state0.x(1,:), GC.cells.num, 1);
state0C.y = repmat(state0.y(1,:), GC.cells.num, 1);
state0C.L = repmat(state0.L(1,:), GC.cells.num, 1);
state0C.K = repmat(state0.K(1,:), GC.cells.num, 1);
state0C.Z_V = repmat(state0.Z_V(1,:), GC.cells.num, 1);
state0C.Z_L = repmat(state0.Z_L(1,:), GC.cells.num, 1);

coarsemodel = getSequential_local(coarsemodel);

W = schedule.control(1).W;
WC = W;
d = pdist2(coarsemodel.transportModel.G.cells.centroids, model.G.cells.centroids([W.cells],:));
[~, c] = min(d);
for wNo = 1:numel(W)
    WC(wNo).cells = c(wNo);
end
coarsemodel.transportNonLinearSolver.LinearSolver = BackslashSolverAD();
coarse = pack(state0C, coarsemodel, 'coarse', 'Coarse');
coarse.SimulatorSetup.schedule.control(1).W = WC;

%% Adaptive reordering

nls = NonLinearSolver();

modeladaptreorder = getSequential_local(model);
modeladaptreorder.pressureModel.extraStateOutput = true;
modeladaptreorder.transportModel = ReorderingModel(modeladaptreorder.transportModel);
modeladaptreorder.transportModel.chunkSize = 5;

modeladaptreorder = AdaptiveSequentialPressureTransportModel(modeladaptreorder.pressureModel, modeladaptreorder.transportModel, GC);

state0.transportState = state0C;
state0.transportState.G = GC;
state0.transportState.pv = modeladaptreorder.transportModel.parent.operators.pv;
state0.transportModel = modeladaptreorder.transportModel;
modeladaptreorder.plotProgress = true;
modeladaptreorder.computeCoarsePressure = true;
modeladaptreorder.refineTol = 1e-2;
modeladaptreorder.coarsenTol = 1e-2;

adaptreorder = pack(state0, modeladaptreorder, 'adapt-reorder', 'Adaptive reordering');

%%

problems = {fim, seq, reorder, adapt, coarse, adaptreorder};

%% Simulate problems

runIx = 6;
for pNo = runIx
    [ok, status] = simulatePackedProblem(problems{pNo});
end

%%

setup = problems{6}.SimulatorSetup;
[ws, st, rep] = simulateScheduleAD(setup.state0, setup.model, setup.schedule);

%% Copyright Notice
%
% <html>
% <p><font size="-1">
% Copyright 2009-2024 SINTEF Digital, Mathematics & Cybernetics.
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
