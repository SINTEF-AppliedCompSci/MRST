mrstModule add dg vem vemmech ad-props ad-core ad-blackoil blackoil-sequential gasinjection mrst-gui reorder matlab_bgl

%%

gravity reset off

n = 10;
l = 1000*meter;
G = computeGeometry(cartGrid([n,1], [1,0.1]*l));
G = computeVEMGeometry(G);
G = computeCellDimensions(G);

rock = makeRock(G, 100*milli*darcy, 1);
fluid = initSimpleADIFluid('phases', 'WO'                   , ...
                           'rho'   , [1000, 1]*kilogram/meter^3, ...
                           'mu'    , [0.5, 0.5]*centi*poise     , ...
                           'n'     , [1, 1]                 );

modelfi = TwoPhaseOilWaterModel(G, rock, fluid);
modelFV = getSequentialModelFromFI(modelfi);
modelDG = modelFV;

%%

time = 2*year;
rate = 2*sum(poreVolume(G, rock))/time;
W = [];
W = addWell(W, G, rock, 1          , 'type', 'rate', 'val', rate    , 'comp_i', [1,0]);
W = addWell(W, G, rock, G.cells.num, 'type', 'bhp' , 'val', 50*barsa, 'comp_i', [1,0]);

if 0
    src = addSource([], 20, rate/3, 'sat', [1,0]);
else
    src = [];
end

dt    = 30*day;
dtvec = rampupTimesteps(time, dt, 0);
schedule = simpleSchedule(dtvec, 'W', W, 'src', src);

sW     = 0.0;
state0 = initResSol(G, 100*barsa, [sW,1-sW]);
% state0.cells = (1:G.cells.num);

%%

degree = [0, 1, 2, 3, 4, 5];
[jt, ot, mt] = deal(Inf);
% 
jt = Inf;
mt = 0.0;
ot = 0.0;
% ot = 0.2;


[wsDG, statesDG, disc] = deal(cell(numel(degree),1));
for dNo = 1:numel(degree)
    disc{dNo} = DGDiscretization(modelDG.transportModel                    , ...
                                    'degree'             , degree(dNo), ...
                                    'basis'              , 'legendre' , ...
                                    'useUnstructCubature', false      , ...
                                    'jumpTolerance'      , jt         , ...
                                    'outTolerance'       , ot         , ...
                                    'meanTolerance'      , mt         , ...
                                    'limiterType'        , 'tvb'      );
    modelDG.transportModel = TransportOilWaterModelDG(G, rock, fluid, ...
                                       'disc'    , disc{dNo}        , ...
                                       'dsMaxAbs', 0.2/(degree(dNo)+1), ...
                                       'nonlinearTolerance', 1e-3);

    state0 = assignDofFromState(modelDG.transportModel.disc, state0);
    [wsDG{dNo}, statesDG{dNo}, rep] = simulateScheduleAD(state0, modelDG, schedule);
    
end

%%

[wsFV, statesFV, rep] = simulateScheduleAD(state0, modelFV, schedule);

%%

close all

for dNo = 1:numel(degree)
    figure
    plotToolbar(G, statesDG{dNo}, 'plot1d', true);
    colormap(jet);
end

figure
plotToolbar(G, statesFV, 'plot1d', true);
colormap(jet)

dsnDG = cellfun(@(d) ['dG(' num2str(d), ')'] , num2cell(degree), 'unif', false);
dsn = horzcat('FV', dsnDG);

plotWellSols({wsFV, wsDG{:}}, schedule.step.val, 'datasetNames', dsn)
% plotWellSols({wsFV, wsDG{:}, wsDGReorder}, schedule.step.val)

%%
close all
figure('position', [-1000, 0, 800, 600])
% h = zeros(numel(degree),1);
saturation = cell(numel(degree), 1);
hold on
for dNo = 1:numel(degree)
    [h(dNo), saturation{dNo}, coords, keep, n] = ...
        plotSaturationDG(disc{dNo}, statesDG{dNo}{1}, 'plot1d', true, 'linew', 2);
end
ll = cellfun(@(d) ['dG(', num2str(d), ')'], num2cell(degree), 'unif', false);
legend(ll);
hold off
axis tight; box on
ylim([-0.1,1.1])

for sNo = 1:numel(statesDG{1})
    for dNo = 1:numel(degree)
        s = saturation{dNo}(statesDG{dNo}{sNo});
        set(h(dNo), 'YData', s);
    end
    pause(0.5);
end

%%

close all

nsteps = 4;
steps = round(linspace(1,numel(schedule.step.val)-10, nsteps));
dNo = 2;

figure('Position', [-2000,0, 2000, 2000]);
nlines = 6;

nclr = 7;

ixDG = 2;
for sNo = 1:nsteps
    
    subplot(2, nsteps, sNo);
    sDG = statesDG{ixDG}{steps(sNo)}.s(:,1);
    sDG = reshape(sDG, [n,n]);
    contourf(sDG, nlines);
    axis equal
    caxis([0,1])
    colormap(summer(nclr))
%     
    subplot(2, nsteps, nsteps + sNo);
    sFV = statesFV{steps(sNo)}.s(:,1);
    sFV = reshape(sFV, [n,n]);
    contourf(sFV, nlines);
    axis equal
    caxis([0,1])
    colormap(summer(nclr))
    
end

% figure('Position', [-2000,0, 1500, 374]);
% for sNo = 1:nsteps
%     sDG = statesDG{ixDG}{steps(sNo)}.s(:,1);
%     sDG = reshape(sDG, [n,n]);
%     sFV = statesFV{steps(sNo)}.s(:,1);
%     sFV = reshape(sFV, [n,n]);
%     
%     subplot(1, nsteps, sNo);
%     yyaxis left
%     contour(sDG, nlines, 'color', 'r');
%     yyaxis right
%     contour(sFV, nlines, 'color', 'b');
%     axis equal
%     
% end


%%

% [jt, ot, mt] = deal(Inf);

[wsDGReorder, statesDGReorder] = deal(cell(numel(degree),1));
for dNo = 1:numel(degree)
    disc = DGDiscretization(modelDG.transportModel, 2, 'degree', degree(dNo), 'basis', 'legendre', 'useUnstructCubature', false, 'jumpTolerance', jt, 'outTolerance', ot, 'meanTolerance', mt);
    modelDG.transportModel = TransportOilWaterModelDG(G, rock, fluid, 'disc', disc);    

    [modelDG.transportModel.extraStateOutput, modelDG.pressureModel.extraStateOutput] = deal(true);
    modelDGReorder = modelDG;
    modelDGReorder.pressureModel.extraStateOutput = true;

    modelDGReorder.transportModel = ReorderingModelDG_ghost(modelDGReorder.transportModel, 'plotProgress', false);

    modelDGReorder.transportModel.chunkSize = 1;
    modelDGReorder.transportModel.parent.extraStateOutput = true;

    
    state0 = assignDofFromState(modelDG.transportModel.disc, state0);
    [wsDGReorder{dNo}, statesDGReorder{dNo}, rep] = simulateScheduleAD(state0, modelDGReorder, schedule);
    
end