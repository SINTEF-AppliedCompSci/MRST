mrstModule add spe10 vem vemmech dg ad-core ad-props ad-blackoil blackoil-sequential gasinjection

%%

[state0, model, schedule]  = setupSPE10_AD('layers', 10);

% Make coarse grid
p = partitionUI(model.G, [20, 40, 1]);
GC = generateCoarseGrid(model.G, p);
GC = coarsenGeometry(GC);
GC = storeInteractionRegionCart(GC);

state0.bfactor = [model.fluid.bW(state0.pressure), model.fluid.bO(state0.pressure)];

% Move wells to center of coarse blocks
W = schedule.control(1).W;
wc = [W.cells];
wc = GC.partition(wc);
[G, map] = refineGrid(GC, GC, model.G, wc);

for wNo = 1:numel(W)
    W(wNo).cells = wc(wNo);
end
scheduleC = schedule;
scheduleC.control(1).W = W;

wx = GC.cells.centroids(wc,:);
for wNo = 1:numel(W)
    d = pdist2(model.G.cells.centroids, wx(wNo,:));
    [~, ix] = min(d);
    W(wNo).cells = ix;
end

schedule.control(1).W = W;

for wNo = 1:numel(W)
    d = pdist2(model.G.cells.centroids, wx(wNo,:));
    [~, ix] = min(d);
    W(wNo).cells = ix;
end

%%

plotGrid(G)

%%

modelSIF = getSequentialModelFromFI(model);

modelASI = AdaptiveSequentialPressureTransportModel(modelSIF.pressureModel, modelSIF.transportModel, G);

if 0
    s = getSmootherFunction('type', 'ilu', 'iterations', 1);

    mssolver = MultiscaleVolumeSolverAD(G, 'getSmoother', s);
    msmodel = MultiscalePressureModel(modelASI.pressureModel.G,...
                                      modelASI.pressureModel.rock,...
                                      modelASI.pressureModel.fluid,...
                                      modelASI.pressureModel, mssolver);
    modelASI.pressureModel = msmodel;
    modelASI.pressureModel.pressureModel.other
    modelASI.pressureModel.pressureModel.extraStateOutput = true;
else
    modelASI.pressureModel.extraStateOutput = true;
end

modelSI = upscaleModelTPFA(model, GC.partition);
modelSI = getSequentialModelFromFI(modelSI);

%%

state0C = upscaleState(modelSI.pressureModel, modelSIF.pressureModel, state0);

%%

if 1
    modelASI.plotProgress = true;
end
[wsASI, statesASI, rep] = simulateScheduleAD(state0, modelASI, schedule);

%%

[wsSI, statesSI, rep] = simulateScheduleAD(state0C, modelSI, scheduleC);

%%

[ws, states, rep] = simulateScheduleAD(state0, model, schedule);

%%

[wsSIF, statesSIF, rep] = simulateScheduleAD(state0, modelSIF, schedule);

%%

plotWellSols({wsSI, wsASI, wsSIF}, cumsum(schedule.step.val), 'datasetNames', {'Coarse', 'Adaptive', 'Reference'});

%%

dt = cumsum(schedule.step.val)/day;
wcut = nan(numel(dt,1),3);
for sNo = 1:numel(statesASI)
    wcut(sNo,:) = [wsSI{sNo}(2).wcut, wsASI{sNo}(2).wcut, wsSIF{sNo}(2).wcut];
end

%%

close all
fig = figure('Position', [-2000, 0, 1500, 600]);
M = struct('cdata',[],'colormap',[]);
d = 10;

clr = lines(3);

for sNo = 1:numel(statesASI)
    
    clf
    subplot(1,3,1)
    hold on
    plotCellData(model.G, statesASI{sNo}.s(:,1), 'edgec', 'none')
    plotGrid(statesASI{sNo}.G, 'facecolor', 'none', 'edgealpha', 0.2);
    hold off
    caxis([0.2,0.8]);
    axis equal tight
    
    subplot(1,3,2)
    plotCellData(modelSI.G, statesSI{sNo}.s(:,1), 'edgec', 'none')
    plotGrid(G, 'facec', 'none', 'edgealpha', 0.2);
    caxis([0.2,0.8]);
    axis equal tight
    
    colormap(jet)
    
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
    legend({'Coarse', 'Adaptive', 'Reference'}, 'location', 'northwest');

    if 0
        rect = [d, d, fig.Position(3:4) - [d,d]];
        M(sNo) = getframe(fig, rect);
    end
    
    pause(0.1)
    
end

%%

pth = mrstPath('dg');
name = 'spe10-refinement';
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

