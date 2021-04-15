mrstModule add dg vem vemmech ad-props ad-core ad-blackoil ...
    blackoil-sequential gasinjection reorder matlab_bgl upr coarsegrid ...
    mrst-gui mrst-utils

mrstVerbose on

%%

n  = 50;
l = 1000;
GF = pebiGrid(l/n, [l,l]);
close all
plotGrid(GF)
axis equal tight

%%

GF = computeGeometry(GF);
GF = computeCellDimensions2(GF);

% Cartesian coarse grid
G_cart = cartGrid([100, 100]);
p_cart = partitionUI(G_cart, [20, 20]);
p_cart = sampleFromBox(GF, reshape(p_cart, G_cart.cartDims));
GC = generateCoarseGrid(GF, p_cart);
GC = coarsenGeometry(GC);
GC = coarsenCellDimensions(GC);
% GC = storeInteractionRegionCart(GC);

%%

close all

mappings = getRefinementMappings(GC, GC, GF, [1, GC.cells.num]);
G        = generateCoarseGrid(GF, mappings.newPartition);
G        = coarsenGeometry(G);
G.cells.refined = mappings.refined;
[G.cells.equal, G.faces.equal] = deal(false);
[GF.cells.equal, GF.faces.equal] = deal(false);


% [G, map1] = refineGrid(GC, GC, GF, [1, GC.cells.num]');
G = coarsenCellDimensions(G);

G.cells.ghost = false(G.cells.num,1);
GF.cells.ghost = false(GF.cells.num,1);


%%

gravity reset off

rock  = makeRock(G, 100*milli*darcy, 1);
rockF = makeRock(GF, 100*milli*darcy, 1);

fluid = initSimpleADIFluid('phases', 'WO'                      , ...
                           'rho'   , [1000, 1]*kilogram/meter^3, ...
                           'mu'    , [0.5, 0.5]*centi*poise    , ...
                           'n'     , [1, 1]                    );

modelFIF = TwoPhaseOilWaterModel(GF, rockF, fluid);
modelSIF = getSequentialModelFromFI(modelFIF);


modelFI  = upscaleModelTPFA(modelFIF, G.partition);
modelSI  = getSequentialModelFromFI(modelFI);
modelASI = AdaptiveSequentialPressureTransportModel(modelSIF.pressureModel, modelSIF.transportModel, G);

modelASI.transportModel.G = coarsenCellDimensions(modelASI.transportModel.G);

[jt, ot, mt] = deal(Inf);
mt = 1e-8;
ot = 1e-8;
degree = 1;
discDG = DGDiscretization(modelASI.transportModel               , ...
                        'degree'               , degree       , ...
                        'basis'                , 'legendre'   , ...
                        'useUnstructCubature'  , true        , ...
                        'jumpTolerance'        , jt           , ...
                        'outTolerance'         , ot           , ...
                        'outLimiter'           , 'orderReduce', ...
                        'meanTolerance'        , mt           , ...
                        'limitAfterConvergence', false        , ...
                        'plotLimiterProgress'  , false        );
transportModelDG = TransportOilWaterModelDG(GF, rockF, fluid, ...
                                   'disc'    , disc        , ...
                                   'dsMaxAbs', 0.2, ...
                                   'nonlinearTolerance', 1e-3);

modelASIDG = AdaptiveSequentialPressureTransportModel(modelSIF.pressureModel, transportModelDG, G);
% modelASIDG.pressureModel = PressureOilWaterModelSemiDG(GF, rockF, fluid);

tm = ReorderingModelDG(transportModelDG);
modelASIDGreorder = AdaptiveSequentialPressureTransportModel(modelASIDG.pressureModel, tm, G);

% modelASIreorder = AdaptiveSequentialPressureTransportModel(modelSIF.pressureModel, transportModel, G);
modelASIDGreorder.transportModel.chunkSize = 50;
modelASIDGreorder.transportModel.parent.extraStateOutput = true;

disc = DGDiscretization(modelSIF.transportModel               , ...
                        'degree'               , degree       , ...
                        'basis'                , 'legendre'   , ...
                        'useUnstructCubature'  , true        , ...
                        'jumpTolerance'        , jt           , ...
                        'outTolerance'         , ot           , ...
                        'outLimiter'           , 'kill', ...
                        'meanTolerance'        , mt           , ...
                        'limitAfterConvergence', false        , ...
                        'plotLimiterProgress'  , false        );
transportModelDG = TransportOilWaterModelDG(GF, rockF, fluid, ...
                                   'disc'    , disc        , ...
                                   'dsMaxAbs', 0.2, ...
                                   'nonlinearTolerance', 1e-3);
modelDG = SequentialPressureTransportModel(modelSIF.pressureModel, transportModelDG);

tm = ReorderingModel(modelSIF.transportModel, 'chunkSize', 50);
G.cells.equal = false;
modelASIreorder = AdaptiveSequentialPressureTransportModel(modelSIF.pressureModel, tm, G);
    

%%

time = 2*year;
rate = 1.5*sum(poreVolume(GF, rockF))/time;

xw = [0,0; l,l];
% xw = G.cells.centroids([1, G.cells.num], :);

W = [];
d = pdist2(G.cells.centroids, xw);

[~, ix] = min(d(:,1));
W = addWell(W, G, rock, ix, 'type', 'rate', 'val', rate    , 'comp_i', [1,0]);
[~, ix] = min(d(:,2));
W = addWell(W, G, rock, ix, 'type', 'bhp' , 'val', 50*barsa, 'comp_i', [1,0]);

dt    = 7*day;
dtvec = rampupTimesteps2(time, dt, 0);
% dtvec = rampupTimesteps(time, dt, 0);
schedule = simpleSchedule(dtvec, 'W', W);

sW     = 0.0;
state0 = initResSol(G, 100*barsa, [sW,1-sW]);
state0.bfactor = [fluid.bW(state0.pressure), fluid.bO(state0.pressure)];

WF = [];
d = pdist2(GF.cells.centroids, xw);

[~, ix] = min(d(:,1));
WF = addWell(WF, GF, rockF, ix, 'type', 'rate', 'val', rate    , 'comp_i', [1,0]);
[~, ix] = min(d(:,2));
WF = addWell(WF, GF, rockF, ix, 'type', 'bhp' , 'val', 50*barsa, 'comp_i', [1,0]);

scheduleF = simpleSchedule(dtvec, 'W', WF);

state0F = initResSol(GF, 100*barsa, [sW,1-sW]);
state0F.bfactor = [fluid.bW(state0F.pressure), fluid.bO(state0F.pressure)];

%%

modelASI.storeGrids                     = true;
modelASI.plotProgress                   = true;
modelASI.pressureModel.extraStateOutput = true;

state0F.transportState        = state0;
state0F.transportState.degree = degree*ones(G.cells.num, 1);
state0F.transportState.G      = G;
state0F.transportState.pv     = modelASI.transportModel.operators.pv;
state0F.transportState        = assignDofFromState(discDG, state0F.transportState);
state0F.transportState        = discDG.updateDofPos(state0F.transportState);
state0F.transportModel        = modelASI.transportModel;

%%

[wsASI, stASI, rep] = simulateScheduleAD(state0F, modelASI, scheduleF);

%%

state0F.transportModel  = modelASIDG.transportModel;
state0F.transportState.G = G;
modelASIDG.plotProgress = true;
[wsASIDG, stASIDG, rep] = simulateScheduleAD(state0F, modelASIDG, scheduleF);

%%

state0F = assignDofFromState(modelDG.transportModel.disc, state0F);
[wsDG, stDG, rep] = simulateScheduleAD(state0F, modelDG, scheduleF);

%%

wsASIreorder = [];

close all
figure('Position', [-1000, 0, 800, 800]);
cmap = mrstColormap();
colormap(cmap);
modelASIDGreorder.transportModel.parent.extraStateOutput = true;
modelASIDGreorder.pressureModel.extraStateOutput = true;
modelASIDGreorder.transportModel.plotProgress = true;
modelASIDGreorder.storeGrids = true;
modelASIDGreorder.plotProgress = false;
modelASIDGreorder.computeCoarsePressure = true;
modelASIDGreorder.transportModel.nonlinearSolver.errorOnFailure = true;
modelASIDGreorder.transportModel.nonlinearSolver.continueOnFailure = true;

state0F.transportModel        = modelASIDGreorder.transportModel;
[wsASIreorder, stASIreorder, rep] = simulateScheduleAD(state0F, modelASIDGreorder, scheduleF);

%%

state0F.transportModel = modelASIreorder.transportModel;
modelASIreorder.plotProgress = true;
modelASIreorder.computeCoarsePressure = true;
[wsASIreorder, stASIreorder, rep] = simulateScheduleAD(state0F, modelASIreorder, scheduleF);

%%

[wsSI, stSI, rep] = simulateScheduleAD(state0, modelSI, schedule);

%%

[wsSIF, stSIF, rep] = simulateScheduleAD(state0F, modelSIF, scheduleF);

%%

close all
cmap = mrstColormap('type', 'wateroil');

%%

ws = {wsSI, wsASIreorder wsASI, wsASIDG, wsSIF};
names    = {'Coarse', ...
            ['Adaptive dG(', num2str(degree), ') reorder'], ...
            'Adaptive dG(0)'                              , ...
            ['Adaptive dG(', num2str(degree), ')']        , ...
            'Fine'};
        
pIx = [1,3:5];
ws = ws(pIx);
names = names(pIx);
        
plotWellSols(ws, scheduleF.step.val, 'datasetNames', names);

%%

figure
plotToolbar(GF, stASI)
axis equal tight
colormap(cmap);

%%

figure
plotToolbar(GF, stASIreorder)
axis equal tight
colormap(cmap);

%%



%%

close all
plotToolbar(G, stDG)
axis equal tight
colormap(jet);

%%

close all
plotToolbar(G, stSI)
axis equal tight
colormap(jet);

%%

dt = cumsum(schedule.step.val)/day;
wcut = nan(numel(dt,1),3);
for sNo = 1:numel(stASI)
    wcut(sNo,:) = [wsSI{sNo}(2).wcut, wsASIDG{sNo}(2).wcut, wsDG{sNo}(2).wcut];
end

%%

close all
fig = figure('Position', [-2000, 0, 1500, 600]);
M = struct('cdata',[],'colormap',[]);
d = 10;

clr = lines(4);
    

% subplot(1,3,3)
% hold on
% wc = wcut;
% wc(2:end,:) = nan;
% xlim([0, dt(end)]);
% ylim([0,1]);
% hold off
% pbaspect([1,1,1])
% box on
% ylabel('Water cut')
% xlabel('Time (days)')
% legend({'Coarse', 'Adaptive', 'Reference'}, 'location', 'northwest');

for sNo = 1:numel(stASI)
    
    subplot(1,3,1)
    cla;
    hold on
    plotCellData(GF, stASI{sNo}.s(:,1), 'edgec', 'none')
    plotGrid(stASI{sNo}.G, 'facecolor', 'none')
    hold off
    caxis([0,1]);
    title('Adaptive dG(1)')
    axis equal tight
    
    subplot(1,3,2)
    plotCellData(G, stSI{sNo}.s(:,1), 'edgec', 'none')
    plotGrid(G, 'facec', 'none', 'edgealpha', 0.2);
    caxis([0,1]);
    title('Coarse dG(0)')
    axis equal tight
    
    colormap(cmap)
    
    subplot(1,3,3)
    wc = wcut;
    wc(sNo+1:end,:) = nan;
    hold on
    plot(dt, wc(:,1), 'linew', 2, 'color', clr(1,:));
    plot(dt, wc(:,2), 'linew', 2, 'color', clr(2,:));
    plot(dt, wc(:,3), '--', 'linew', 4, 'color', clr(3,:));
    hold off
    pbaspect([1,1,1])
    xlim([0, dt(end)]);
    ylim([0,1]);
    box on
    ylabel('Water cut')
    xlabel('Time (days)')
    legend({'Coarse dG(0)', 'Adaptive dG(1)' 'Fine dG(1)'}, 'location', 'northwest');
%     for wNo = 1:3
%         h(wNo).YData = wc(:,wNo);
%     end

    if 1
        rect = [d, d, fig.Position(3:4) - [d,d]];
        M(sNo) = getframe(fig, rect);
    end
    
    pause(0.05)
    
end

%%

pth = mrstPath('dg');
name = 'refinement-2';
duration = 10;
vo = VideoWriter(fullfile(pth, name));
vo.FrameRate = numel(stASI)/duration;
open(vo);

writeVideo(vo, M);

close(vo)

%%

rock = makeRock(GF, 1,1);
T = computeTrans(GF, rock);
p = partitionMETIS(GF, T, 50);
GC = generateCoarseGrid(GF, p);
GC = coarsenGeometry(GC);
GC = addCoarseCenterPoints(GC);
GC = coarsenCellDimensions(GC);

%%

close all
cNo = 50;
% plotGrid(G, cNo, 'facec', [1,1,1]*0.3);
p = (GC.partition == cNo);
plotGrid(GF, p, 'facec', [1,1,1]*0.8);

hold on
dx = GC.cells.dx(cNo,:)/2;
xc = GC.cells.centroids(cNo,:);
xMin = GC.cells.xMin(cNo,:);
x = [xc(1) - dx(1), xc(2) - dx(2);
     xc(1) + dx(1), xc(2) - dx(2); 
     xc(1) + dx(1), xc(2) + dx(2); 
     xc(1) - dx(1), xc(2) + dx(2); 
     xc(1) - dx(1), xc(2) - dx(2)];
plot(x(:,1), x(:,2), 'k', 'linew', 2);
plot(xc(1), xc(2), '.k', 'markerSize', 15)

f = 1.2;
axis([xc(1) - dx(1)*f, xc(1) + dx(1)*f, xc(2) - dx(2)*f, xc(2) + dx(2)*f]);

axis equal

%% Copyright Notice
%
% <html>
% <p><font size="-1">
% Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.
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
