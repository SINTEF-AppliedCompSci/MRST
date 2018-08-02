mrstModule add dg vem vemmech ad-props ad-core ad-blackoil blackoil-sequential gasinjection mrst-gui

%%

gravity reset off
% gravity reset on; gravity([0,-9.81]);



n = 20;
l = 1000;
G = computeGeometry(cartGrid([n,n], [l,l]*meter));
G = computeVEMGeometry(G);
G = computeCellDimensions(G);

rock = makeRock(G, 100*milli*darcy, 1);
fluid = initSimpleADIFluid('phases', 'WO'                   , ...
                           'rho'   , [1000, 1]*kilogram/meter^3, ...
                           'mu'    , [0.5, 0.5]*centi*poise     , ...
                           'n'     , [1, 1]                 );
% fluid.krW = @(s) fluid.krW(s).*(s>=0 & s<=1) + fluid.krW(1).*(s>1);
% fluid.krO = @(s) fluid.krO(s).*(s>=0 & s<=1) + fluid.krO(1).*(s>1);

modelfi = TwoPhaseOilWaterModel(G, rock, fluid);
modelFV = getSequentialModelFromFI(modelfi);
modelDG = modelFV;

%%

time = 2*year;
rate = 1.2*sum(poreVolume(G, rock))/time;
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

degree = [0, 1];
% degree = 1;

[jt, ot, mt] = deal(Inf);
% 
% jt = 0.2;
% ot = 0.1;
% mt = 0;

[wsDG, statesDG] = deal(cell(numel(degree),1));
for dNo = 1:numel(degree)
    disc = DGDiscretization(modelDG.transportModel, 2, 'degree', degree(dNo), 'basis', 'legendre', 'useUnstructCubature', false, 'jumpTolerance', jt, 'outTolerance', ot, 'meanTolerance', mt);
    modelDG.transportModel = TransportOilWaterModelDG(G, rock, fluid, 'disc', disc);    

    state0 = assignDofFromState(modelDG.transportModel.disc, state0);
    [wsDG{dNo}, statesDG{dNo}, rep] = simulateScheduleAD(state0, modelDG, schedule);
end

%%

[modelDG.transportModel.extraStateOutput, modelDG.pressureModel.extraStateOutput] = deal(true);
modelDGreorder = modelDG;
modelDGreorder.pressureModel.extraStateOutput = true;

modelDGreorder.transportModel = ReorderingModelDG_ghost(modelDGreorder.transportModel, 'plotProgress', true);

modelDGreorder.transportModel.chunkSize = 6;
modelDGreorder.transportModel.parent.extraStateOutput = true;


[wsDGReorder, statesDGReorder, repDGReorder] = simulateScheduleAD(state0, modelDGreorder, schedule);


%%

[wsFV, statesFV, rep] = simulateScheduleAD(state0, modelFV, schedule);

%%

close all

nclr = 6;
for dNo = 1:numel(degree)
    figure
    plotToolbar(G, statesDG{dNo});
    colormap(jet);
end

figure
plotToolbar(G, statesFV);
colormap(jet)

figure
plotToolbar(G, statesDGReorder);
colormap(jet)

plotWellSols({wsFV, wsDG{:}, wsDGReorder}, schedule.step.val)

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
