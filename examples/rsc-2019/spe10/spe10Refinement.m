mrstModule add spe10 vem vemmech dg ad-core ad-props ad-blackoil ...
    blackoil-sequential gasinjection incomp msrsb matlab_bgl coarsegrid ...
    mrst-gui reorder weno

mrstVerbose on

%% Common stuff

dz = -3;
pw = @(G,W) plot3(G.cells.centroids([W.cells], 1), ...
                  G.cells.centroids([W.cells], 2), ...
                  G.cells.centroids([W.cells], 3) + dz, ...
                 'ok', 'markerSize', 8, 'markerFaceColor', 'w', 'lineWidth', 2);

baseName = 'spe10';
% location = 'home';
location = 'work';
switch location
    case 'work'
        dataDir  = fullfile('/media/strene/806AB4786AB46C92/mrst-dg/rsc-2019', baseName);
    case 'home'
        dataDir  = fullfile('/home/oysteskl/Documents/phd/repositories/mrst-dg/simulation-output/rsc-2019', baseName);
end
            
%% Set up base model

[state0, model, schedule]  = setupSPE10_AD('layers', 70);
G        = model.G;
rock     = model.rock;
fluid     = model.fluid;
% schedule = schedule;
schedule.step.val = rampupTimesteps(sum(schedule.step.val), 20*day, 8);
schedule.step.control = ones(numel(schedule.step.val), 1);
G        = computeCellDimensions2(G);
G.equal  = true;

%% Construct coarse grid

T = computeTrans(G, rock);
p = partitionMETIS(G, T, 800);
GC = generateCoarseGrid(G, p);
GC = coarsenGeometry(GC);
GC = coarsenCellDimensions(GC, 'useQhull', true);
GC = storeInteractionRegionCart(GC);

%% Make new wells

time = 2000*day;
bhp  = 50*barsa;
rate = 0.25*sum(poreVolume(G, rock))/time;

WF = [];
WF = verticalWell(WF, G, rock, 1, 1, []  , ...
                                   'type'  , 'rate', ...
                                   'val'   , rate  , ...
                                   'comp_i', [1,0] , ...
                                   'name'  , 'Inj' );
WF = verticalWell(WF, G, rock, 42, 220, [], ...
                                   'type'  , 'bhp' , ...
                                   'val'   , bhp   , ...
                                   'comp_i', [1,0] , ...
                                   'name'  , 'Prod');
schedule.control(1).W = WF;

%% Make coarse grid with refinement around wells

xw = G.cells.centroids([WF.cells], :);

d      = pdist2(GC.cells.centroids, xw);
refine = any(d < 50*meter,2);

mappings = getRefinementMappings(GC, GC, G, refine);
GC        = generateCoarseGrid(G, mappings.newPartition);
GC.cells.refined = mappings.refined;
GC        = coarsenGeometry(GC);
GC        = coarsenCellDimensions(GC, 'useQhull', true);
GC.equal  = false;

% %% Make well struct for the coarse grid
% 
% wc = GC.partition([WF.cells]);
W = WF;
for wNo = 1:numel(W)
    W(wNo).cells = GC.partition(W(wNo).cells);
end
% schedule = schedule;
% schedule.control(1).W  = W;

%% Plot setup

close all
plotToolbar(G, rock);
plotGrid(GC, 'facec', 'none', 'edgec', [1,1,1]*0.25, 'linew', 1.5);
hold on
pw(GC, W);
colormap(pink)
axis equal tight; ax = gca;

%% Set up models

degree = [0,1];
jt = 0.1; ot = 1e-6; mt = 1e-6;

% Fully implicit fine model
G.cells.equal = true;
G.faces.equal = false;
GC.cells.equal = false;
GC.faces.equal = false;
modelFI = TwoPhaseOilWaterModel(G, rock, model.fluid);

[modelDG, modelDGreorder, modelDGadapt, modelDGadaptReorder] = deal(cell(numel(degree),1));

for dNo = 1:numel(degree)
    
    % Assign sequential models
    [modelDG{dNo}, modelDGreorder{dNo}, modelDGadapt{dNo}, modelDGadaptReorder{dNo}] ...
        = deal(getSequentialModelFromFI(modelFI));

    % Set up discretization
    disc = DGDiscretization(modelDG{dNo}.transportModel, ...
                   'degree'               , degree(dNo), ...
                   'basis'                , 'legendre' , ...
                   'useUnstructCubature'  , true       , ...
                   'jumpTolerance'        , jt         , ...
                   'outTolerance'         , ot         , ...
                   'outLimiter'           , 'kill'     , ...
                   'meanTolerance'        , mt         );
    % Transport model
    tmodel = TransportOilWaterModelDG(G, rock, model.fluid, ...
                                'disc'              , disc, ...
                                'dsMaxAbs'          , 0.1 , ...
                                'nonlinearTolerance', 1e-3);
                            
    % Reordering model
    tmodelReorder = ReorderingModelDG(tmodel, 'chunkSize', 1);
    
    % Regular model
    modelDG{dNo}.transportModel = tmodel;
    
    % Reordering model
    modelDGreorder{dNo}.transportModel = tmodelReorder;
    
    % Adaptive model
    modelDGadapt{dNo} ...
        = AdaptiveSequentialPressureTransportModel(modelDG{dNo}.pressureModel, tmodel, GC);
                    
    % Adaptive reordering mdoel
    modelDGadaptReorder{dNo} ...
        = AdaptiveSequentialPressureTransportModel(modelDG{dNo}.pressureModel, tmodelReorder, GC);
    
end

modelRef = getSequentialModelFromFI(modelFI);

%%

state0.bfactor = [fluid.bW(state0.pressure), fluid.bO(state0.pressure)];
state0C = initResSol(GC, state0.pressure(1), state0.s(1,:));
state0C.bfactor = [fluid.bW(state0C.pressure), fluid.bO(state0C.pressure)];

%% Set up problems

names = {'base', 'reorder', 'adapt', 'adapt-reorder'};
problems = cell(numel(degree), numel(modelDG{1}));
for dNo = 1:numel(degree)
    
    tmodel = modelDG{dNo}.transportModel;
    state0         = assignDofFromState(tmodel.disc, state0);
    state0.bfactor = [fluid.bW(state0.pressure), fluid.bO(state0.pressure)];
    
    tmodel = modelDGadapt{dNo}.transportModel;
    state0C = assignDofFromState(tmodel.disc, state0C);
    pv = tmodel.operators.pv;
    d = fullfile(dataDir, ['dg', num2str(degree(dNo))]);
    mdls = {modelDG{dNo}, modelDGreorder{dNo}, modelDGadapt{dNo}, modelDGadaptReorder{dNo}};
    
    for mNo = 1:numel(mdls)
        
        if isa(mdls{mNo}, 'AdaptiveSequentialPressureTransportModel')
            tmodel  = mdls{mNo}.transportModel;
            state0.transportState = state0C;
            state0.transportState.pv = pv;
            state0.transportState.G = GC;
            state0.transportModel = tmodel;
            mdls{mNo}.computeCoarsePressure = true;
            mdls{mNo}.plotProgress = true;
        end
        mdls{dNo}.pressureModel.extraStateOutput = true;
        problems{dNo,mNo} = packSimulationProblem(state0, mdls{mNo}, schedule, baseName, 'Directory', d, 'name', names{mNo});
    end
end

%%

mdlIx = 1;
for dNo = 1:2%numel(degree)
    for mNo = mdlIx
        [ok, status] = simulatePackedProblem(problems{dNo,mNo});
    end
end

%%

setup = problems{1,4}.SimulatorSetup;
setup.model.plotProgress = true;
[ws, st, rep] = simulateScheduleAD(setup.state0, setup.model, setup.schedule);

%%



models = {modelSIF, modelSI    , modelASI  , modelASIReorder   , modelDGF, modelDG    , modelASIDG, modelASIDGReorder};
names  = {'fv'    , 'fv-coarse', 'fv-adapt', 'fv-adapt-reorder', 'dg'    , 'dg-coarse', 'dg-adapt', 'dg-adapt-reorder'};
problems = cell(size(models));

state0         = initResSol(G, state0.pressure(1), state0.s(1,:));
state0         = assignDofFromState(disc, state0);
state0.bfactor = [fluid.bW(state0.pressure), fluid.bO(state0.pressure)];
state0.bfactor = [fluid.bW(state0.pressure), fluid.bO(state0.pressure)];
state0         = assignDofFromState(disc, state0);
state0C        = initResSol(GC, state0.pressure(1), state0.s(1,:));


for mNo = 1:numel(models)

    st0 = state0;
    sch = schedule;
    mdl = models{mNo};
    mdl.extraStateOutput = true;
    if isa(mdl, 'AdaptiveSequentialPressureTransportModel')
        st0.transportState        = state0;
        st0.transportState.degree = degree*ones(GC.cells.num, 1);
        st0.transportState.G      = GC;
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

runIx = 6;
for pNo = runIx
    [ok, status] = simulatePackedProblem(problems{pNo});
end

%%

plotIx = 6;
for pNo = plotIx
   [ws, st, rep] = getPackedSimulatorOutput(problems{pNo});
end

%%

close all
figure
plotToolbar(G, st)
axis equal tight
colormap(pink)

%%

%%

problems{1} = packSimulationProblem(state0, modelSIF, schedule);

%%

% Coarse state
state0 = initResSol(GC, state0.pressure(1), state0.s(1,:));




%%
state0.transportState        = state0;
state0.transportState.degree = degree*ones(GC.cells.num, 1);
state0.transportState.G      = GC;
state0.transportState.pv     = modelASI.transportModel.operators.pv;

state0.transportModel = modelASI.transportModel;
modelASI.plotProgress  = true;
close all
[wsASI, stASI, rep] = simulateScheduleAD(state0, modelASI, schedule);

%%

state0.transportModel = modelASI.transportModel;
[wsSIF, stSIF, rep] = simulateScheduleAD(state0, modelSIF, schedule);

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
state0 = initResSol(GC, state0.pressure(1), state0.s(1,:));
state0 = assignDofFromState(disc, state0);
[wsSIDG, stSIDG, rep] = simulateScheduleAD(state0, modelDGF, subschedule);

%%

modelASIDG.storeGrids                     = true;
modelASIDG.plotProgress                   = true;
modelASIDG.pressureModel.extraStateOutput = true;
modelASIDG.transportModel.dsMaxAbs = 0.05;
modelASIDG.transportModel.disc.limitAfterConvergence = true;
% modelASIDG.reupdatePressure = true;

state0.transportState        = state0;
state0.transportState.degree = degree*ones(GC.cells.num, 1);
state0.transportState.G      = GC;
state0.transportState.pv     = modelASIDG.transportModel.operators.pv;
state0.transportState        = disc.assignDofFromState(state0.transportState);
state0.transportState        = disc.updateDofPos(state0.transportState);
state0.transportModel        = modelASIDG.transportModel;



nls = NonLinearSolver('useLineSearch', false);
modelASIDG.transportNonLinearSolver = nls;

[wsASIDG, stASIDG, rep] = simulateScheduleAD(state0, modelASIDG, schedule);

%%

state0 = assignDofFromState(disc, state0);
[wsSIFDG, stSIFDG, rep] = simulateScheduleAD(state0, modelSIFDG, schedule);

%%

modelASI.storeGrids                     = true;
modelASI.plotProgress                   = true;
modelASI.pressureModel.extraStateOutput = true;

state0 = initResSol(GC, state0.pressure(1), state0.s(1,:));
state0.transportState        = state0;
state0.transportState.degree = degree*ones(GC.cells.num, 1);
state0.transportState.G      = GC;
state0.transportState.pv     = modelASIDG.transportModel.operators.pv;
state0.transportState        = disc.assignDofFromState(state0.transportState);
state0.transportState        = disc.updateDofPos(state0.transportState);
state0.transportModel        = modelASI.transportModel;

[wsASI, statesASI, rep] = simulateScheduleAD(state0, modelASI, schedule);

%%

[ws, states, rep] = simulateScheduleAD(state0, model, schedule);


%%

plotWellSols({wsSI, wsASI, wsSIF}, cumsum(schedule.step.val), 'datasetNames', {'Coarse', 'Adaptive', 'Fine'});

%%

close all;
pos = [-1000, 0, 600, 800];

%%

figure('position', pos)
plotToolbar(G, stASI);
cmap = mrstColormap();
colormap(cmap);
axis equal tight; box on;

%%

figure('position', pos)
plotToolbar(G, stSIF);
cmap = mrstColormap();
colormap(cmap);
axis equal tight; box on;

%%

figure('position', pos)
plotToolbar(GC, stSI);
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
    h(mNo) = plotCellData(GC, states{mNo}{1}.s(:,1), 'edgec', 'none');
    colorbar('Location', 'southoutside');
    axis equal off
    text(x(mNo), 20, titles{mNo}, 'fontsize', 25, 'color', 'w'); 
    caxis([0,1])
end

subplot(1,3,1);
h(3) = plotCellData(GC, states{1}{1}.degree, 'edgec', 'none');
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
plotFaces(G, faces, 'facealpha', 0.2);
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

%%

clr = lines(3);
close all
figure; hold on
fNo   = floor(rand*(GC.faces.num-1)) + 1;
faces = GC.faces.fconn(GC.faces.connPos(fNo):GC.faces.connPos(fNo+1)-1);
plotFaces(G, faces, 'facealpha', 0.2);
axis equal tight

xf = GC.faces.centroids(fNo,:);
vec1 = GC.faces.coordSys{1}(fNo,:);
vec2 = GC.faces.coordSys{2}(fNo,:);
n    = GC.faces.normals(fNo,:)./GC.faces.areas(fNo);
mag = 5;
quiver3(xf(:,1), xf(:,2), xf(:,3), vec1(1), vec1(2), vec1(3), mag, 'linew', 2, 'color', clr(1,:));
quiver3(xf(:,1), xf(:,2), xf(:,3), vec2(1), vec2(2), vec2(3), mag, 'linew', 2, 'color', clr(1,:));
quiver3(xf(:,1), xf(:,2), xf(:,3), n(1)   , n(2)   , n(3)   , mag , 'linew', 2, 'color', clr(2,:));

[~, x] = disc.surfaceCubature.getCubature(fNo, 'face');

% n = G.faces.nodes(mcolon(G.faces.nodePos(fNo), G.faces.nodePos(fNo+1)-1));
% x = G.nodes.coords(n,:);
plot3(x(:,1), x(:,2), x(:,3), '.', 'markerSize', 15);

% e = G.faces.edges(G.faces.edgePos(fNo):G.faces.edgePos(fNo+1)-1);
% n = G.edges.nodes(mcolon(G.edges.nodePos(e), G.edges.nodePos(e+1)-1));
% x = G.nodes.coords(n,:);
% plot3(x(:,1), x(:,2), x(:,3), 'o', 'markerSize', 15);
% 
x = GC.nodes.coords(GC.faces.nodes(GC.faces.nodePos(fNo):GC.faces.nodePos(fNo+1)-1),:);
xx = [x; x(1,:)];
patch(xx(:,1), xx(:,2), xx(:,3), 'r', 'facealpha', 0.2);
view(3)

%%

clr = lines(3);
close all
figure; hold on
fNo   = floor(rand*(GC.faces.num-1)) + 1;
faces = GC.faces.fconn(GC.faces.connPos(fNo):GC.faces.connPos(fNo+1)-1);
plotFaces(G, faces, 'facealpha', 0.2);
axis equal tight

xf = GC.faces.centroids(fNo,:);
vec1 = GC.faces.coordSys{1}(fNo,:);
vec2 = GC.faces.coordSys{2}(fNo,:);
n    = GC.faces.normals(fNo,:)./GC.faces.areas(fNo);
mag = 5;
quiver3(xf(:,1), xf(:,2), xf(:,3), vec1(1), vec1(2), vec1(3), mag, 'linew', 2, 'color', clr(1,:));
quiver3(xf(:,1), xf(:,2), xf(:,3), vec2(1), vec2(2), vec2(3), mag, 'linew', 2, 'color', clr(1,:));
quiver3(xf(:,1), xf(:,2), xf(:,3), n(1)   , n(2)   , n(3)   , mag , 'linew', 2, 'color', clr(2,:));

disc.surfaceCubature
% 
[~, x] = disc.surfaceCubature.getCubature(fNo, 'face');
% 
% n = G.faces.nodes(mcolon(G.faces.nodePos(fNo), G.faces.nodePos(fNo+1)-1));
% x = G.nodes.coords(n,:);
plot3(x(:,1), x(:,2), x(:,3), '.', 'markerSize', 15);

% e = G.faces.edges(G.faces.edgePos(fNo):G.faces.edgePos(fNo+1)-1);
% n = G.edges.nodes(mcolon(G.edges.nodePos(e), G.edges.nodePos(e+1)-1));
% x = G.nodes.coords(n,:);
% plot3(x(:,1), x(:,2), x(:,3), 'o', 'markerSize', 15);
% 
x = GC.nodes.coords(GC.faces.nodes(GC.faces.nodePos(fNo):GC.faces.nodePos(fNo+1)-1),:);
xx = [x; x(1,:)];
patch(xx(:,1), xx(:,2), xx(:,3), 'r', 'facealpha', 0.2);
view(3)