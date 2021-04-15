mrstModule add dg vem vemmech ad-props ad-core ad-blackoil upr...
    blackoil-sequential mrst-gui reorder matlab_bgl ad-eor trust-region ...
    weno deckformat
mrstVerbose on
gravity reset off

%%

[G, rock, fluid, state0, schedule] = setupSaigupBC('layers', 1:3);
[G.cells.equal, G.faces.equal] = deal(false);

%%

baseName = 'saigup';
dataDir  = fullfile(mrstPath('dg'), 'examples', 'siamgs-19', 'data', baseName);
pack = @(name, state0, model, schedule) ...
            packSimulationProblem(state0, model, schedule, baseName, ...
                                               'Directory', dataDir, ...
                                               'Name'     , name   );

[mt, ot] = deal(1e-6);
%%

modelFI   = TwoPhaseOilWaterModel(G, rock, fluid);
problemFI = pack('fi', state0, modelFI, schedule);

%%

simulatePackedProblem(problemFI);

%%

modelSI   = getSequentialModelFromFI(modelFI);
modelSI.transportModel.nonlinearTolerance = 1e-7;
problemSI = pack('si', state0, modelSI, schedule);

%%

simulatePackedProblem(problemSI);

%%

modelDG = getSequentialModelFromFI(modelFI);
modelDG.transportModel ...
    = TransportOilWaterModelDG(G, rock, fluid, ...
                                'dsMaxAbs'     , 0.2, ...
                                'degree'       , 0  , ...
                                'meanTolerance', mt , ...
                                'outTolerance' , ot , ...
                                'jumpTolerance', Inf, ...
                                'useUnstructCubature', false);
modelDG.transportModel.nonlinearTolerance = 1e-3;
% modelDG.transportModel.disc.degree0 = ones(G.cells.num,1);
% modelDG.transportModel.disc.degree0(c) = 0;
state0     = assignDofFromState(modelDG.transportModel.disc, state0);
problemDG0 = pack('dg-0', state0, modelDG, schedule);

%%

simulatePackedProblem(problemDG0);

%%

modelDG = getSequentialModelFromFI(modelFI);
modelDG.transportModel ...
    = TransportOilWaterModelDG(G, rock, fluid, ...
                                'dsMaxAbs'     , 0.2, ...
                                'degree'       , 1  , ...
                                'meanTolerance', mt , ...
                                'outTolerance' , ot , ...
                                'jumpTolerance', Inf, ...
                                'useUnstructCubature', false);
modelDG.transportModel.nonlinearTolerance = 1e-3;
% modelDG.transportModel.disc.degree0 = ones(G.cells.num,1);
% modelDG.transportModel.disc.degree0(c) = 0;
[state0.posFlux, state0.negFlux] = deal(false);
state0     = assignDofFromState(modelDG.transportModel.disc, state0);

problemDG1 = pack('dg-1', state0, modelDG, schedule);

%%

simulatePackedProblem(problemDG1);

%%

modelWENO = HigherOrderOilWaterModel(G, rock, fluid);
weno      = WENODiscretization(modelWENO, G.griddim, ...
                             'includeBoundary'     , true, ...
                             'interpolateReference', true);
modelWENO.FluxDiscretization.saturationDiscretization = weno;
modelWENO.FluxDiscretization.relPermDiscretization = weno;
modelWENO.FluxDiscretization.discritizeRelPerm     = false;
modelWENO.FluxDiscretization.discritizeViscosity   = false;
% modelWENO.dpMaxRel = 0.1;
problemWENO = pack('weno', state0, modelWENO, schedule);

%%

simulatePackedProblem(problemWENO);

%%

% setup = problemDG0.SimulatorSetup;
setup = problemDG1Reorder.SimulatorSetup;
% setup = problemSI.SimulatorSetup;
[ws, st, rep] = simulateScheduleAD(setup.state0, setup.model, setup.schedule);

%%

[~, stDG0 , repDG0 ] = getPackedSimulatorOutput(problemDG0);
[~, stDG1 , repDG1 ] = getPackedSimulatorOutput(problemDG1);
[~, stWENO, repWENO] = getPackedSimulatorOutput(problemWENO);

%%

close all

pba  = [1,3,0.4];
azel = [-100, 35];
pos  = [0,0,1400,600];
lpos = [1500, 9500, 0];
% lpos = [1,1,2];

cmap = winter();
cmap = cmap(end:-1:1,:);
sMin = 0.25;

figDir = fullfile(mrstPath('dg'), 'examples', 'siamgs-19', 'fig', 'saigup');
% savepng = @(name) print(fullfile(figDir, name), '-dpng', '-r400');
savepng = @(name) []

%%

close all

hdg0  = figure('name', 'dG(0)', 'position', pos);
hdg1  = figure('name', 'dG(1)', 'position', pos);
hweno = figure('name', 'WENO', 'position', pos);

% stepNo = [10:20:120];
stepNo = [70, 100, 120, 180];

for sNo = stepNo
    
    if problemDG0.OutputHandlers.states.numelData >= sNo
        figure(hdg0); clf;
        s = stDG0{sNo}.s(:,1);
        plotCellData(G, s, s>sMin);
        plotGrid(G, 'facec', 'none', 'edgealpha', 0.2);
        axis tight
        pbaspect(pba);
        view(azel)
        colormap(cmap);
        light('position', lpos, 'style', 'local');
        caxis([sMin,1]);
        ax = gca; [ax.XTick, ax.YTick, ax.ZTick] = deal([]);
        drawnow();pause(0.1);
        savepng(['saigup-sat-dg0-', num2str(sNo)]);
    end

    figure(hdg1); clf;
    if problemDG1.OutputHandlers.states.numelData >= sNo
        s = stDG1{sNo}.s(:,1);
        plotCellData(G, s, s>sMin);
        plotGrid(G, 'facec', 'none', 'edgealpha', 0.2);
        axis tight
        pbaspect(pba);
        view(azel)
        colormap(cmap);
        caxis([sMin,1]);
        light('position', lpos, 'style', 'local');
        ax = gca; [ax.XTick, ax.YTick, ax.ZTick] = deal([]);
        drawnow();pause(0.1);
        savepng(['saigup-sat-dg1-', num2str(sNo)]);
    end

    if problemWENO.OutputHandlers.states.numelData >= sNo
        figure(hweno); clf;
        s = stDG1{sNo}.s(:,1);
        plotCellData(G, s, s>sMin);
        plotGrid(G, 'facec', 'none', 'edgealpha', 0.2);
        axis tight
        pbaspect(pba);
        view(azel)
        colormap(cmap);
        caxis([sMin,1]);
        light('position', lpos, 'style', 'local');
        ax = gca; [ax.XTick, ax.YTick, ax.ZTick] = deal([]);
        drawnow();pause(0.1);
        savepng(['saigup-sat-weno-', num2str(sNo)]);
    end
end

%%

close all
[Gfull, rockfull] = setupSaigupBC();

%%

figure('name', 'perm', 'position', pos);
perm = rockfull.perm(:,1);
plotCellData(Gfull, log10(perm), 'edgealpha', 0.2);
pbaspect(pba);
view(azel)
[hb, hh] = colorbarHist(perm, [min(perm), max(perm)], 'South', 500, true);
colormap(pink);
d = 10;
hb.Position = hb.Position + [0,0,0,0.02];
hb.TickLabels = round((10.^hb.Ticks)/(milli*darcy)*100)/100;
hb.FontSize = 12;

savepng('saigup-perm');

%%

close all

figure('name', 'grid', 'position', pos);
[nf, cNo] = max(diff(Gfull.cells.facePos));
f = Gfull.cells.faces(Gfull.cells.facePos(cNo):Gfull.cells.facePos(cNo+1)-1);
n = Gfull.faces.neighbors(f,:);
n = n(:);
plotGrid(Gfull, 'facec', 'none', 'edgealpha', 0.2)
plotGrid(Gfull, n, 'edgealpha', 0.2);
axis tight
pbaspect(pba);

%%

close all

figure('name', 'grid', 'position', [0,0,800,400]);

% plotGrid(Gfull, cNo);
plotGrid(Gfull, n, 'edgealpha', 0.3, 'facec', 'none');

axis tight
axis equal off

view([-130,30])
% light('position', [1800,1500, 1000], 'style', 'local');
light('position', [1500,1800, 2000], 'style', 'local');
light('position', [1500,1800, 2000], 'style', 'local');
gr = [1,1,1]*0.8;

ffmin = min(Gfull.faces.areas(f));
ffmax = max(Gfull.faces.areas(f));
clr = lines(numel(f));
for fNo = 1:numel(f)
    ff = f(fNo);
    if Gfull.faces.areas(ff) < 2000
        cc = clr(fNo,:);
    else
        cc = gr;
    end
    plotFaces(Gfull, ff, 'facecolor', cc)
end

savepng('saigup-cells');

%%

close all

figure('name', 'grid', 'position', [0,0,800,400]);

% plotGrid(Gfull, cNo);
% plotGrid(Gfull, n, 'edgealpha', 0.3, 'facec', 'none');

axis tight
axis equal off

view([-47,0])
% light('position', [1800,1500, 1000], 'style', 'local');
light('position', [1500,1800, 2000], 'style', 'local');
light('position', [1500,1800, 2000], 'style', 'local');
gr = [1,1,1]*0.8;

ffmin = min(Gfull.faces.areas(f));
ffmax = max(Gfull.faces.areas(f));
clr = lines(numel(f));
for fNo = 1:numel(f)
    ff = f(fNo);
    if Gfull.faces.areas(ff) < 2000
        cc = clr(fNo,:);
    else
        cc = gr;
    end
    plotFaces(Gfull, ff, 'facecolor', cc)
end

xc = Gfull.cells.centroids(cNo,:);
hold on
plot3(xc(1), xc(2), xc(3), '.k', 'markerSize', 20)

savepng('saigup-centroid');

%%

close all

cc = randi(Gfull.cells.num,1)
hold on
plotGrid(Gfull, cc, 'facec', 'none');
xc = Gfull.cells.centroids(cc,:);
plot3(xc(1), xc(2), xc(3), '.k', 'markerSize', 20)
axis equal tight

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
