mrstModule add spe10 vem vemmech dg ad-core ad-props ad-blackoil ...
    blackoil-sequential gasinjection incomp msrsb matlab_bgl coarsegrid ...
    mrst-gui reorder

mrstVerbose on

%% Common stuff

dz = -3;
pw = @(G,W) plot3(G.cells.centroids([W.cells], 1), ...
                  G.cells.centroids([W.cells], 2), ...
                  G.cells.centroids([W.cells], 3) + dz, ...
                 'ok', 'markerSize', 8, 'markerFaceColor', 'w', 'lineWidth', 2);

baseName = 'spe10';
dataDir  = fullfile('/media/strene/806AB4786AB46C92/mrst-dg/rsc-2019', baseName);

%% Set up base model

[state0F, model, schedule]  = setupSPE10_AD('layers', 70);
GF        = model.G;
rockF     = model.rock;
fluid     = model.fluid;
scheduleF = schedule;
scheduleF.step.val = rampupTimesteps2(sum(schedule.step.val), 20*day, 8);
scheduleF.step.control = ones(numel(scheduleF.step.val), 1);
GF        = computeCellDimensions2(GF);
GF.equal  = true;

%% Construct coarse grid

T = computeTrans(GF, rockF);
p = partitionMETIS(GF, T, 800);
GC = generateCoarseGrid(GF, p);
GC = coarsenGeometry(GC);
GC = coarsenCellDimensions(GC, 'useQhull', true);
GC = storeInteractionRegionCart(GC);

%% Make new wells

time = 2000*day;
bhp  = 50*barsa;
rate = 0.25*sum(poreVolume(GF, rockF))/time;

WF = [];
WF = verticalWell(WF, GF, rockF, 1, 1, []  , ...
                                   'type'  , 'rate', ...
                                   'val'   , rate  , ...
                                   'comp_i', [1,0] , ...
                                   'name'  , 'Inj' );
WF = verticalWell(WF, GF, rockF, 42, 220, [], ...
                                   'type'  , 'bhp' , ...
                                   'val'   , bhp   , ...
                                   'comp_i', [1,0] , ...
                                   'name'  , 'Prod');
scheduleF.control(1).W = WF;

%% Make coarse grid with refinement around wells

xw = GF.cells.centroids([WF.cells], :);

d      = pdist2(GC.cells.centroids, xw);
refine = any(d < 50*meter,2);

mappings = getRefinementMappings(GC, GC, GF, refine);
G        = generateCoarseGrid(GF, mappings.newPartition);
G.cells.refined = mappings.refined;
G        = coarsenGeometry(G);
G        = coarsenCellDimensions(G, 'useQhull', true);
G.equal  = false;

%% Make well struct for the coarse grid

wc = GC.partition([WF.cells]);
W = WF;
for wNo = 1:numel(W)
    W(wNo).cells = G.partition(W(wNo).cells);
end
schedule = scheduleF;
schedule.control(1).W  = W;

%% Plot setup

close all
plotToolbar(GF, rockF);
plotGrid(G, 'facec', 'none', 'edgec', [1,1,1]*0.25, 'linew', 1.5);
hold on
pw(G, W);
colormap(pink)
axis equal tight; ax = gca;

%% Set up models

degree = 1;
jt = 0.1; ot = 1e-6; mt = 1e-6;

% Fully implicit fine model
modelFIF = TwoPhaseOilWaterModel(GF, rockF, model.fluid);
% Sequential fine model
modelSIF = getSequentialModelFromFI(modelFIF);                                             
% Fully implicit coarse model
modelFI = upscaleModelTPFA(modelFIF, G.partition);
% Sequential coarse model
modelSI = getSequentialModelFromFI(modelFI);
% Adaptive sequential model
modelASI = AdaptiveSequentialPressureTransportModel(modelSIF.pressureModel, modelSIF.transportModel, G);
% Adaptive reordering FV model
transportModelReorder = ReorderingModel(modelSIF.transportModel);
transportModelReorder.chunkSize = 100;
modelASIReorder = AdaptiveSequentialPressureTransportModel(modelSIF.pressureModel, transportModelReorder, G);

% Fine DG model
modelDGF = getSequentialModelFromFI(modelFIF);
discDGF = DGDiscretization(modelDGF.transportModel         , ...
                        'degree'               , degree    , ...
                        'basis'                , 'legendre', ...
                        'useUnstructCubature'  , true      , ...
                        'jumpTolerance'        , jt        , ...
                        'outTolerance'         , ot        , ...
                        'outLimiter'           , 'kill'    , ...
                        'meanTolerance'        , mt        );
transportModelDGF = TransportOilWaterModelDG(GF, rockF, model.fluid, ...
                                      'disc'              , discDGF, ...
                                      'dsMaxAbs'          , 0.1    , ...
                                      'nonlinearTolerance', 1e-3   );
modelDGF.transportModel = transportModelDGF;
% Coarse DG model
modelDG = getSequentialModelFromFI(modelFI);
modelDG.transportModel.G = coarsenCellDimensions(modelDG.transportModel.G);
discDG = DGDiscretization(modelDG.transportModel          , ...
                       'degree'               , degree    , ...
                       'basis'                , 'legendre', ...
                       'useUnstructCubature'  , true      , ...
                       'jumpTolerance'        , jt        , ...
                       'outTolerance'         , ot        , ...
                       'outLimiter'           , 'kill'    , ...
                       'meanTolerance'        , mt        );
transportModelDG = TransportOilWaterModelDG(GF, rockF, model.fluid, ...
                                     'disc'              , discDG , ...
                                     'dsMaxAbs'          , 0.1    , ...
                                     'nonlinearTolerance', 1e-3   );
modelDG.transportModel = transportModelDG;
% Adaptive DG model
modelASIDG = AdaptiveSequentialPressureTransportModel(modelDGF.pressureModel, modelDGF.transportModel, G);
% Adaptive reordering DG model
transportModelDGReorder = ReorderingModelDG(modelDGF.transportModel, 'chunkSize', 10);
modelASIDGReorder = AdaptiveSequentialPressureTransportModel(modelDGF.pressureModel, transportModelDGReorder, G);

%% Set up problems

models = {modelSIF, modelSI    , modelASI  , modelASIReorder   , modelDGF, modelDG    , modelASIDG, modelASIDGReorder};
names  = {'fv'    , 'fv-coarse', 'fv-adapt', 'fv-adapt-reorder', 'dg'    , 'dg-coarse', 'dg-adapt', 'dg-adapt-reorder'};
problems = cell(size(models));

state0          = initResSol(G, state0F.pressure(1), state0F.s(1,:));
state0.bfactor  = [fluid.bW(state0.pressure), fluid.bO(state0.pressure)];
state0F.bfactor = [fluid.bW(state0F.pressure), fluid.bO(state0F.pressure)];
state0F         = assignDofFromState(discDGF, state0F);


for mNo = 1:numel(models)

    st0 = state0F;
    sch = scheduleF;
    mdl = models{mNo};
    mdl.extraStateOutput = true;
    if isa(mdl, 'AdaptiveSequentialPressureTransportModel')
        st0.transportState        = state0;
        st0.transportState.degree = degree*ones(G.cells.num, 1);
        st0.transportState.G      = G;
        tm = mdl.transportModel;
        if mdl.isReordering
            tm = tm.parent;
            mdl.computeCoarsePressure = true;
        end
        st0.transportState.pv     = tm.operators.pv;
        st0.transportModel        = mdl.transportModel;
    elseif isCoarseGrid(mdl.transportModel.G)
        st0 = state0;
        sch = schedule;
    end
    problems{mNo} = packSimulationProblem(st0, mdl, sch, baseName, ...
                                          'Directory', dataDir   , ...
                                          'name'     , names{mNo});
end

%% Simulate problems

runIx = [1,2,3,5,6,7];
for pNo = runIx
    [ok, status] = simulatePackedProblem(problems{pNo});
end

%%

setup = problems{5}.SimulatorSetup;

[ws, st, rep] = simulateScheduleAD(setup.state0, setup.model, setup.schedule);

%%

problems{1} = packSimulationProblem(state0F, modelSIF, scheduleF);

%%

% Coarse state
state0 = initResSol(G, state0F.pressure(1), state0F.s(1,:));




%%
state0F.transportState        = state0;
state0F.transportState.degree = degree*ones(G.cells.num, 1);
state0F.transportState.G      = G;
state0F.transportState.pv     = modelASI.transportModel.operators.pv;

state0F.transportModel = modelASI.transportModel;
modelASI.plotProgress  = true;
close all
[wsASI, stASI, rep] = simulateScheduleAD(state0F, modelASI, scheduleF);

%%

state0F.transportModel = modelASI.transportModel;
[wsSIF, stSIF, rep] = simulateScheduleAD(state0F, modelSIF, scheduleF);

%%

[wsSI, stSI, rep] = simulateScheduleAD(state0, modelSI, schedule);

%%

ix = 1:75;
subschedule = schedule;
subschedule.step.val = subschedule.step.val(ix);
subschedule.step.control = subschedule.step.control(ix);

modelDGF.transportModel.dsMaxAbs = 0.1;
modelDGF.transportModel.disc.meanTolerance = mt;
modelDGF.transportModel.disc.outTolerance = ot;
modelDGF.transportModel.disc.jumpTolerance = jt;

% state0.degree = degree*ones(G.cells.num,1);
state0 = initResSol(G, state0F.pressure(1), state0F.s(1,:));
state0 = assignDofFromState(discDG, state0);
[wsSIDG, stSIDG, rep] = simulateScheduleAD(state0, modelDGF, subschedule);

%%

modelASIDG.storeGrids                     = true;
modelASIDG.plotProgress                   = true;
modelASIDG.pressureModel.extraStateOutput = true;
modelASIDG.transportModel.dsMaxAbs = 0.05;
modelASIDG.transportModel.disc.limitAfterConvergence = true;
% modelASIDG.reupdatePressure = true;

state0F.transportState        = state0;
state0F.transportState.degree = degree*ones(G.cells.num, 1);
state0F.transportState.G      = G;
state0F.transportState.pv     = modelASIDG.transportModel.operators.pv;
state0F.transportState        = discDG.assignDofFromState(state0F.transportState);
state0F.transportState        = discDG.updateDofPos(state0F.transportState);
state0F.transportModel        = modelASIDG.transportModel;



nls = NonLinearSolver('useLineSearch', false);
modelASIDG.transportNonLinearSolver = nls;

[wsASIDG, stASIDG, rep] = simulateScheduleAD(state0F, modelASIDG, scheduleF);

%%

state0F = assignDofFromState(discDGF, state0F);
[wsSIFDG, stSIFDG, rep] = simulateScheduleAD(state0F, modelSIFDG, scheduleF);

%%

modelASI.storeGrids                     = true;
modelASI.plotProgress                   = true;
modelASI.pressureModel.extraStateOutput = true;

state0 = initResSol(G, state0F.pressure(1), state0F.s(1,:));
state0F.transportState        = state0;
state0F.transportState.degree = degree*ones(G.cells.num, 1);
state0F.transportState.G      = G;
state0F.transportState.pv     = modelASIDG.transportModel.operators.pv;
state0F.transportState        = discDG.assignDofFromState(state0F.transportState);
state0F.transportState        = discDG.updateDofPos(state0F.transportState);
state0F.transportModel        = modelASI.transportModel;

[wsASI, statesASI, rep] = simulateScheduleAD(state0F, modelASI, scheduleF);

%%

[ws, states, rep] = simulateScheduleAD(state0, model, schedule);


%%

plotWellSols({wsSI, wsASI, wsSIF}, cumsum(scheduleF.step.val), 'datasetNames', {'Coarse', 'Adaptive', 'Fine'});

%%

close all;
pos = [-1000, 0, 600, 800];

%%

figure('position', pos)
plotToolbar(GF, stASI);
cmap = mrstColormap();
colormap(cmap);
axis equal tight; box on;

%%

figure('position', pos)
plotToolbar(GF, stSIF);
cmap = mrstColormap();
colormap(cmap);
axis equal tight; box on;

%%

figure('position', pos)
plotToolbar(G, stSI);
cmap = mrstColormap();
colormap(cmap);
axis equal tight; box on;

%%

dt = cumsum(schedule.step.val)/day;
wcut = nan(numel(dt,1),3);
for sNo = 1:numel(stASI)
    wcut(sNo,:) = [wsSI{sNo}(2).wcut, wsASI{sNo}(2).wcut, wsSIF{sNo}(2).wcut];
end

%%

close all
fig = figure('Position', [-2000, 0, 1500, 600]);
M = struct('cdata',[],'colormap',[]);
d = 10;

clr = lines(3);

for sNo = 1:numel(stASI)
    
    clf
    subplot(1,3,1)
    hold on
    plotCellData(modelASI.G, stASI{sNo}.s(:,1), 'edgec', 'none')
    plotGrid(stASI{sNo}.G, 'facecolor', 'none', 'edgealpha', 0.2);
    hold off
    caxis([0.2,0.8]);
    axis equal tight
    title('Adaptive')
    
    subplot(1,3,2)
    plotCellData(modelSIF.G, stSIF{sNo}.s(:,1), 'edgec', 'none')
    plotGrid(modelSIF.G, 'facecolor', 'none', 'edgealpha', 0.2);
    caxis([0.2,0.8]);
    axis equal tight
    title('Fine')
    
    colormap(cmap)
    
    subplot(1,3,3)
    hold on
    wc = wcut;
    wc(sNo+1:end,:) = nan;
    plot(dt, wc(:,1), 'linew', 2, 'color', clr(1,:));
    plot(dt, wc(:,2), 'linew', 2, 'color', clr(2,:));
    plot(dt, wc(:,3), '--', 'linew', 4, 'color', clr(3,:));
    xlim([0, dt(end)]);
    ylim([0,1]);
    hold off
    pbaspect([1,1,1])
    box on
    ylabel('Water cut')
    xlabel('Time (days)')
    legend({'Coarse', 'Adaptive', 'Fine'}, 'location', 'northwest');

    if 1
        rect = [d, d, fig.Position(3:4) - [d,d]];
        M(sNo) = getframe(fig, rect);
    end
    
    pause(0.1)
    
end

%%

pth = mrstPath('dg');
name = 'spe10-refinement-2';
duration = 10;
vo = VideoWriter(fullfile(pth, name));
vo.FrameRate = numel(stASI)/duration;
open(vo);

writeVideo(vo, M);

close(vo)

%%

close all

figure('position', [-2000, 0, 400, 800]);
plotToolbar(model.G, statesASI);

figure('position', [-2000, 0, 400, 800]);
plotToolbar(model.G, statesSIF);


%%

close all

fig = figure('Position', [-2000, 0, 1500, 1000]);

states = {statesDG, statesFV};
titles = {'dG(1)', 'FV'};

x = [250, 295];
for mNo = 1:numel(states)
    subplot(1,3,mNo+1)
    h(mNo) = plotCellData(G, states{mNo}{1}.s(:,1), 'edgec', 'none');
    colorbar('Location', 'southoutside');
    axis equal off
    text(x(mNo), 20, titles{mNo}, 'fontsize', 25, 'color', 'w'); 
    caxis([0,1])
end

subplot(1,3,1);
h(3) = plotCellData(G, states{1}{1}.degree, 'edgec', 'none');
colormap jet
caxis([0,1]);
colorbar('Location', 'southoutside');
axis equal off
text(150, 20, 'dG degree', 'fontsize', 25, 'color', 'w'); 

% set(fig, 'Units', 'pixels');
% pos = get(fig, 'Position');
% set(fig, 'Units', 'normalized');
M = struct('cdata',[],'colormap',[]);

for sNo = 1:numel(schedule.step.val)
    
    for mNo = 1:numel(states)
        h(mNo).CData = states{mNo}{sNo}.s(:,1);
    end
    
    h(3).CData = states{1}{sNo}.degree;
    
    pause(0.01);
    drawnow
%     dx = 10;
%     dy = 10;
    rect = [0, 0, fig.Position(3:4)];
%     rect = pos;
    M(sNo) = getframe(fig, rect);
    
end

%%

pth = mrstPath('dg');
name = 'spe10';
duration = 10;
vo = VideoWriter(fullfile(pth, name));
vo.FrameRate = numel(states{1})/duration;
open(vo);

writeVideo(vo, M);

close(vo)

%%

clr = lines(3);
close all
figure; hold on
fNo   = floor(rand*(GC.faces.num-1)) + 1;
faces = GC.faces.fconn(GC.faces.connPos(fNo):GC.faces.connPos(fNo+1)-1);
plotFaces(GF, faces, 'facealpha', 0.2);
axis equal tight

xf = GC.faces.centroids(fNo,:);
vec1 = GC.faces.coordSys{1}(fNo,:);
vec2 = GC.faces.coordSys{2}(fNo,:);
n    = GC.faces.normals(fNo,:)./GC.faces.areas(fNo);
quiver3(xf(:,1), xf(:,2), xf(:,3), vec1(1), vec1(2), vec1(3), 20, 'linew', 2, 'color', clr(1,:));
quiver3(xf(:,1), xf(:,2), xf(:,3), vec2(1), vec2(2), vec2(3), 20, 'linew', 2, 'color', clr(1,:));
quiver3(xf(:,1), xf(:,2), xf(:,3), n(1)   , n(2)   , n(3)   , 20 , 'linew', 2, 'color', clr(2,:));


n = GC.faces.nodes(mcolon(GC.faces.nodePos(fNo), GC.faces.nodePos(fNo+1)-1));
x = GC.nodes.coords(n,:);
plot3(x(:,1), x(:,2), x(:,3), '.', 'markerSize', 15);

e = GC.faces.edges(GC.faces.edgePos(fNo):GC.faces.edgePos(fNo+1)-1);
n = GC.edges.nodes(mcolon(GC.edges.nodePos(e), GC.edges.nodePos(e+1)-1));
x = GC.nodes.coords(n,:);
plot3(x(:,1), x(:,2), x(:,3), 'o', 'markerSize', 15);

x = GC.nodes.coords(GC.faces.nodes(GC.faces.nodePos(fNo):GC.faces.nodePos(fNo+1)-1),:);
xx = [x; x(1,:)];
patch(xx(:,1), xx(:,2), xx(:,3), 'r', 'facealpha', 0.2);
view(3)


