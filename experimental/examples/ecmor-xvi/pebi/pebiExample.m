mrstModule add dg vem vemmech ad-props ad-core ad-blackoil blackoil-sequential gasinjection reorder matlab_bgl upr mrst-gui spe10 gasinjection
mrstVerbose on

%%

if 1

    xmax = [500, 1000, 500]*meter;
%     xmax = [500, 1000, 500]/2*meter;
%     xmax = [100, 300, 100]*meter;
%     n = [5, 10, 5];
    n = [10, 20, 10];

    xx = cell(3,1);
    d = 5*meter;
    for dNo = 1:3
        xx{dNo} = linspace(d,xmax(dNo)-d,n(dNo));
    end

    [x,y,z] = ndgrid(xx{:});
    x = [x(:), y(:), z(:)];
    npts = size(x,1);
    
    rng(1);
    d = xmax./n*0.2;
    x = x + randn(npts,3).*d;

    bnd = [0      , 0      , 0      ;
           xmax(1), 0      , 0      ;
           xmax(1), xmax(2), 0      ;
           0      , xmax(2), 0      ;
           0      , 0      , xmax(3);
           xmax(1), 0      , xmax(3);
           xmax(1), xmax(2), xmax(3);
           0      , xmax(2), xmax(3)];
    bnd = delaunayTriangulation(bnd);
    G = clippedPebi3D(x, bnd);

    G = computeVEMGeometry(G);
    G = computeCellDimensions(G);

else
    
    load('grid.mat')
    
end
%% 

dy = 115;
dx = 10;
rock  = makeRock(G, 100*milli*darcy, 0.4);

if 1
    [~, m, ~] = setupSPE10_AD('layers', 1:20);
    perm = reshape(m.rock.perm(:,1), [m.G.cartDims(1:2), 20]);
    % perm = perm((1:2*n(1)) + dx,(1:2*n(2)) + dy,:);
    rock.perm = sampleFromBox(G, perm);
%     rock.perm = repmat(rock.perm, 1, 3);
%     rock.perm(:,1) = rock.perm(:,1)*10;

    poro = reshape(m.rock.poro(:,1), [m.G.cartDims(1:2), 20]);
    % poro = poro((1:2*n(1)) + dx,(1:2*n(2)) + dy,:);
    rock.poro = sampleFromBox(G, poro);
end
%%

fluid = initSimpleADIFluid('phases', 'WO', ...
                           'n'     , [1,1], ...
                           'mu'    , [0.5,0.5]*centi*poise, ...
                           'rho'   , [1,1]);

time = 4*year;
xw = [xmax(1)/2, 0, xmax(3)/2; xmax(1)/2, xmax(2), xmax(3)/2];
isInj = [true false];
rate = 1.5*sum(poreVolume(G, rock))/time;
bhp = 50*barsa;
type = {'bhp', 'bhp'};
val = [800*barsa, bhp];

W = [];
for wNo = 1:size(xw,1)
    
    d = sqrt(sum((G.cells.centroids - xw(wNo,:)).^2, 2));
    c = find(d == min(d)); c = c(1);
    W = addWell(W, G, rock, c, 'type', type{wNo}, 'val', val(wNo), 'comp_i', [1,0]);
    
end

dt = 60*day;
dtvec = rampupTimesteps(time, dt, 0);
schedule = simpleSchedule(dtvec, 'W', W);

%%

close all
plotToolbar(G, rock);
% plotGrid(G);
plotWell(G, W);
axis equal tight
view(3)

%%

dataDir = '/media/strene/806AB4786AB46C92/mrst-dg/pebi';
getOutHandler = @(name) ResultHandler('dataDirectory', dataDir, ...
                                   'dataFolder'   , name   , ...
                                   'dataPrefix'   , 'state' , ...
                                   'cleardir'     , false  );
                               
getRepHandler = @(name) ResultHandler('dataDirectory', dataDir, ...
                                   'dataFolder'   , name   , ...
                                   'dataPrefix'   , 'rep', ...
                                   'cleardir'     , false  );
                               
%%

modelFI = TwoPhaseOilWaterModel(G, rock, fluid);
modelFV = getSequentialModelFromFI(modelFI);
sW     = 0.0;
state0 = initResSol(G, 100*barsa, [sW, 1-sW]);

%%

[wsFV, statesFV, repFV] = simulateScheduleAD(state0, modelFV, schedule);

%%


modelDG = modelFV;
[modelDG.transportModel.extraStateOutput, ...
 modelDG.pressureModel.extraStateOutput] = deal(true);

[jt, ot, mt] = deal(Inf);
ot = 1e-1;
% jt = 0.3;
% ot = 0;


degree = [0,1];
% [wsDG, statesDG, repDG, wsDGReorder, statesDGReorder, repDGReorder] = deal(cell(numel(degree),1));
runIx = 2;
for dNo = runIx
    
    disc   = DGDiscretization(modelDG.transportModel, G.griddim, ...
                              'degree'             , degree(dNo), ...
                              'basis'              , 'legendre' , ...
                              'useUnstructCubature', true       , ...
                              'jumpTolerance'      , jt         , ...
                              'outTolerance'       , ot         , ...
                              'meanTolerance'      , mt         );
    modelDG.transportModel = TransportOilWaterModelDG(G, rock, fluid, ...
                                'disc'              , disc, ...
                                'nonlinearTolerance', 1e-3, ...
                                'dsMaxAbs'          , 0.1);    
    state0 = assignDofFromState(modelDG.transportModel.disc, state0);
% 
    oh = getOutHandler(['dg',num2str(degree(dNo)),'-2']);
    rh = getRepHandler(['dg',num2str(degree(dNo)),'-2']);
%     [oh, rh] = deal([]);
%     [wsDG{dNo}, statesDG{dNo}, repDG{dNo}] ...
%         = simulateScheduleAD(state0, modelDG, schedule, ...
%                              'OutputHandler', oh, ...
%                              'ReportHandler', rh);
%     
    modelDGReorder = modelDG;
    modelDGReorder.transportModel ...
        = ReorderingModelDG_ghost(modelDGReorder.transportModel, ...
            'plotProgress'      , false, ...
            'plotAfterTimestep' , true , ...
            'chunkSize'         , 300  , ...
            'nonlinearTolerance', 1e-3 );
        
    modelDGReorder.transportModel.parent.extraStateOutput = true;
    modelDGReorder.transportModel.parent.nonlinearTolerance = 1e-3;
    
    oh = getOutHandler(['dg',num2str(degree(dNo)), '-reorder-2']);
    rh = getRepHandler(['dg',num2str(degree(dNo)), '-reorder-2']);
    [oh, rh] = deal([]);
    [wsDGReorder{dNo}, statesDGReorder{dNo}, repDGReorder{dNo}] ...
        = simulateScheduleAD(state0, modelDGReorder, schedule, ...
                             'OutputHandler', oh, ...
                             'ReportHandler', rh);
    
end

%%

pth = fullfile(mrstPath('dg'), 'examples', 'ecmor-xvi', 'pebi', 'fig');

if 0
    savepng = @(name) print(fullfile(pth, name), '-dpng', '-r300');
    saveeps = @(name) print(fullfile(pth, name), '-depsc');
else
    savepng = @(name) [];
    saveeps = @(name) [];
end

%%

close all

stateNo = [5,20,50];
% stateNo = [1:6];
pos = [0,0,800,400];
azel = [120, 15];
gr = [1,1,1]*0.5;

for dNo = 1:numel(degree)
    for sNo = stateNo
        
        figure('Position', pos, 'name', ['dG(', num2str(degree(dNo)), ')'])
        it = repDGReorder{dNo}.ControlstepReports{sNo}.StepReports{1}.NonlinearReport{1}.TransportSolver.StepReports{1}.NonlinearReport{1}.Iterations;
        ii = it > 0;
        plotGrid(G, 'facec', 'none');
%         plotCellData(G, it(it>0), it>0);
        plotGrid(G, it>0, 'facecolor', gr);
        plotWell(G, W, 'color', 'black', 'height', 170)
        axis equal off
        view(azel)
        savepng(['pebi-solved-dg', num2str(degree(dNo)), '-', num2str(sNo)]);
        
        figure('position', pos, 'name', ['dG(', num2str(degree(dNo)), ')']);
        s = statesDGReorder{dNo}{sNo}.s(:,1);
        ii = s > 0.5;
        plotCellData(G, s(ii), ii);
        plotWell(G, W, 'color', 'black', 'height', 170)
        caxis([0,1]);
        axis equal off
        view(azel)
        savepng(['pebi-sat-dg', num2str(degree(dNo)), '-', num2str(sNo)]);
        
    end
end

%%

close all

wNo = 2;
dtt = cumsum(dtvec)/year;
figure('Position', pos)
hold on
for dNo = 1:numel(degree)
    wcut = cellfun(@(w) w(wNo).wcut, wsDG{dNo});
    plot(dtt, wcut, 'linewidth', 2)
end
axis([0 dtt(end) 0 1])
legend({'dG(0)', 'dG(1)'}, 'location', 'northwest');
box on
xlabel('Time (years)');
ylabel('Watercut');
ax = gca;
ax.FontSize = 15;
saveeps('pebi-wcut')

%%

close all

figure('pos', pos)
plotCellData(G, log10(rock.perm));
plotWell(G, W, 'color', 'black', 'height', 170)
logColorbar('location', 'southoutside')
view(azel)
axis equal off
savepng('pebi-perm');

figure('pos', pos)
plotCellData(G, rock.poro);
plotWell(G, W, 'color', 'black', 'height', 170)
colorbar('location', 'southoutside')
view(azel)
axis equal off
savepng('pebi-poro');


%%

plotCellData(G, statesDGReorder{1}.order)
hold on
x = G.cells.centroids(:,1:2);


close all
%%

close all
for dNo = 1:numel(degree)
    figure('name', ['dG(', num2str(degree(dNo)), ')']);
    plotToolbar(G, statesDG{dNo});
    axis equal tight
end

% for dNo = 1:numel(degree)
%     figure('name', ['Reordered dG(', num2str(degree(dNo)), ')']);
%     plotToolbar(G, statesDGReorder{dNo});
%     axis equal tight
% end


plotWellSols({wsDG{:}, wsDGReorder{:}});
% plotWellSols(wsDG);

%%
figure; hold on
colors = jet(G.cells.num);
for cNo = 1:15:G.cells.num
    s = cellfun(@(x) x.s(cNo, 1), statesDG{1});
    s2 = cellfun(@(x) x.s(cNo, 1), statesDG{2});
    
    plot(s, '--', 'color', colors(cNo, :));
    plot(s2, '-', 'color', colors(cNo, :));
end


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
