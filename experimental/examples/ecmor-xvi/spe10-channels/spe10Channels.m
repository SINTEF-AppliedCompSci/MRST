mrstModule add dg vem vemmech ad-props ad-core ad-blackoil blackoil-sequential gasinjection reorder matlab_bgl upr mrst-gui spe10
mrstVerbose on

%%

[state0, modelFI, schedule] = setupSPE10_AD('layers', 50);
G = modelFI.G;
G = computeVEMGeometry(G);
G = computeCellDimensions(G);
modelFI.G = G;
rock = modelFI.rock;
fluid = modelFI.fluid;
close all

time = 4*year;
rate = 0.2*sum(poreVolume(G, rock))/time;
bhp = 275*barsa;

W = [];
W = verticalWell(W, G, rock, 31, 1, [], ...
                 'type', 'rate', ...
                 'val', rate, ...
                 'comp_i', [1,0]);
W = verticalWell(W, G, rock, 35, 220, [], ...
                 'type', 'bhp', ...
                 'val', bhp, ...
                 'comp_i', [1,0]);
             
dt = 20*day;
dtvec = rampupTimesteps(time, dt);
schedule = simpleSchedule(dtvec, 'W', W);

%%

plotToolbar(modelFI.G, modelFI.rock);
plotWell(G, W);
axis equal tight

%%

modelFV = getSequentialModelFromFI(modelFI);
modelDG = modelFV;
[modelDG.transportModel.extraStateOutput, ...
 modelDG.pressureModel.extraStateOutput  ] = deal(true);

%%

[wsFV, statesFV, repFV] = simulateScheduleAD(state0, modelFV, schedule);

%%

dataDir = '/media/strene/806AB4786AB46C92/mrst-dg/spe10-channels';
getOutHandler = @(name) ResultHandler('dataDirectory', dataDir, ...
                                   'dataFolder'   , name   , ...
                                   'dataPrefix'   , 'state' , ...
                                   'cleardir'     , false  );
                               
getRepHandler = @(name) ResultHandler('dataDirectory', dataDir, ...
                                   'dataFolder'   , name   , ...
                                   'dataPrefix'   , 'rep', ...
                                   'cleardir'     , false  );
                               
%%

modelWENO = HigherOrderOilWaterModel(G, rock, fluid);
stencil =  'local';

weno = WENODiscretization(modelWENO, G.griddim, 'includeBoundary', true,...
    'stencilType', stencil, 'interpolateReference', true);
modelWENO.FluxDiscretization.saturationDiscretization = weno;
modelWENO.FluxDiscretization.relPermDiscretization = weno;
modelWENO.FluxDiscretization.discritizeRelPerm = false;
modelWENO.FluxDiscretization.discritizeViscosity = false;

ohWENO = getOutHandler('weno');
rhWENO = getRepHandler('weno');
[wsWENO, statesWENO, repWENO] = simulateScheduleAD(state0, modelWENO, schedule, ...
                                    'ReportHandler', rhWENO, ...
                                    'OutputHandler', ohWENO);

%%

[jt, ot, mt] = deal(Inf);
ot = 1e-3;

degree = [0,1];
% dsMaxAbs = 0.05;
dsMaxAbs = 0.1;
[wsDG, statesDG, repDG, wsDGReorder, statesDGReorder, repDGReorder] = deal(cell(numel(degree),1));

runIx = 1:2;
for dNo = runIx
    
    disc   = DGDiscretization(modelDG.transportModel, G.griddim, ...
                             'degree'             , degree(dNo), ...
                             'basis'              , 'legendre' , ...
                             'useUnstructCubature', true       ,  ...
                             'jumpTolerance'      , jt         , ...
                             'outTolerance'       , ot         , ...
                             'meanTolerance'      , mt         );
    modelDG.transportModel = TransportOilWaterModelDG(G, rock, fluid, ...
                                'disc'              , disc, ...
                                'nonlinearTolerance', 1e-3, ...
                                'dsMaxAbs'          , dsMaxAbs);
    state0 = assignDofFromState(modelDG.transportModel.disc, state0);
    
    ohDG = getOutHandler(['dg', num2str(degree(dNo))]);
    rhDG = getRepHandler(['dg', num2str(degree(dNo))]);
%     [wsDG{dNo}, statesDG{dNo}, repDG{dNo}] ...
%         = simulateScheduleAD(state0, modelDG, schedule, ...
%                              'OutputHandler', ohDG, ...
%                              'Reporthandler', rhDG);
                         
    modelDGReorder = modelDG;
    modelDGReorder.transportModel ...
        = ReorderingModelDG_ghost(modelDGReorder.transportModel, ...
                                  'plotProgress'      , false, ...
                                  'plotAfterTimestep' , true , ...
                                  'chunkSize'         , 1000 , ...
                                  'nonlinearTolerance', 1e-3 );
    modelDGReorder.transportModel.parent.nonlinearTolerance = 1e-3;
    modelDGReorder.transportModel.parent.extraStateOutput = true;
    
    ohDG = getOutHandler(['dg', num2str(degree(dNo)), '-reorder']);
    rhDG = getRepHandler(['dg', num2str(degree(dNo)), '-reorder']);
    [wsDGReorder{dNo}, statesDGReorder{dNo}, repDGReorder{dNo}] ...
        = simulateScheduleAD(state0, modelDGReorder, schedule, ...
                             'OutputHandler', ohDG, ...
                             'Reporthandler', rhDG);
%     
end

%%

pth = fullfile(mrstPath('dg'), 'examples', 'ecmor-xvi', 'pebi', 'fig');

if 1
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
azel = [6,44];
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
        plotCellData(G, statesDGReorder{dNo}{sNo}.s(:,1));
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


plotWellSols({wsFV, wsDG{2}});
% plotWellSols(wsDG);


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
