mrstModule add dg vem vemmech ad-props ad-core ad-blackoil ...
    blackoil-sequential vista mrst-gui reorder matlab_bgl weno
mrstVerbose on;

%%

gravity reset off
% gravity reset on; gravity([0,-9.81]);



n = 10;
l = 1000;
G = computeGeometry(cartGrid([n,n], [l,l]*meter));
G = computeVEMGeometry(G);
G = computeCellDimensions(G);
[G.cells.equal, G.faces.equal] = deal(true);

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

[bc, src, W] = deal([]);
if 0
    src = addSource(src, 20, rate/3, 'sat', [1,0]);
end
if 1
    [f, c] = boundaryFaces(G);
    f = f(c == 1);
    bc = addBC(bc, f, 'flux', rate/2, 'sat', [1,0]);
else
    W = addWell(W, G, rock, 1, 'type', 'rate' , 'val', rate, 'comp_i', [1,0]);
end
W = addWell(W, G, rock, G.cells.num, 'type', 'bhp' , 'val', 50*barsa, 'comp_i', [1,0]);

dt    = 30*day;
dtvec = rampupTimesteps(time, dt, 0);
schedule = simpleSchedule(dtvec, 'W', W, 'src', src, 'bc', bc);

sW     = 0.0;
state0 = initResSol(G, 100*barsa, [sW,1-sW]);
% state0.cells = (1:G.cells.num);

%%

degree = 0:5;


% degree = [4];
degree = [0,1,2];
[jt, ot, mt] = deal(Inf);
% 
% jt = Inf;
mt = 1e-6;
ot = 1e-3;
% ot = 0.2;


[wsDG, statesDG, disc] = deal(cell(numel(degree),1));
nls = NonLinearSolver('maxIterations', 25, 'useLinesearch', false);
for dNo = 1:numel(degree)
    disc{dNo} = DGDiscretization(modelDG.transportModel                   , ...
                                    'degree'               , degree(dNo)  , ...
                                    'basis'                , 'legendre'   , ...
                                    'useUnstructCubature'  , true        , ...
                                    'jumpTolerance'        , jt           , ...
                                    'outTolerance'         , ot           , ...
                                    'outLimiter'           , 'orderReduce', ...
                                    'meanTolerance'        , mt           , ...
                                    'limitAfterConvergence', false         , ...
                                    'plotLimiterProgress'  , false        );
    modelDG.transportModel = TransportOilWaterModelDG(G, rock, fluid, ...
                                       'disc'    , disc{dNo}        , ...
                                       'dsMaxAbs', 0.2, ...
                                       'nonlinearTolerance', 1e-3);
    
    modelDG.transportModel.conserveOil   = true;
    modelDG.transportModel.conserveWater = false;
    
    modelDG.pressureModel = PressureOilWaterModelSemiDG(G, rock, fluid, ...
                                       'disc'    , disc{dNo}        );
    modelDG.transportNonLinearSolver = nls;

    state0 = assignDofFromState(modelDG.transportModel.disc, state0);
    [wsDG{dNo}, statesDG{dNo}, rep] = simulateScheduleAD(state0, modelDG, schedule);
end

%%

[wsFV, statesFV, rep] = simulateScheduleAD(state0, modelFV, schedule);

%%

close all

dsnDG = cellfun(@(d) ['dG(' num2str(d), ')'] , num2cell(degree), 'unif', false);
dsn = horzcat('FV', dsnDG);

for dNo = 1:numel(degree)
    figure('name', dsnDG{dNo})
    plotToolbar(G, statesDG{dNo});
    colormap(jet);
end

figure
plotToolbar(G, statesFV);
colormap(jet)

plotWellSols({wsFV, wsDG{:}}, schedule.step.val, 'datasetNames', dsn)
% plotWellSols({wsFV, wsDG{:}, wsDGReorder}, schedule.step.val)

%%
close all
figure('position', [-1000, 0, 800, 600])

azel = [107, 16];
pba = [5,5,1];
dNo = 3;
[h, saturation, coords, keep, n] = plotSaturationDG(disc{dNo}, statesDG{dNo}{1}, 'edgecolor', 'none');
view(azel);
pbaspect(pba)
colormap(jet)
light('position', [1100, 500, 2], 'style', 'local');
axis tight; box on

st = statesDG{dNo};
for sNo = 1:numel(st)
    s = saturation(st{sNo});
    s(~keep) = nan;
    s = reshape(s', [n, n]);
    h.ZData = s;
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

%% Copyright Notice
%
% <html>
% <p><font size="-1">
% Copyright 2009-2023 SINTEF Digital, Mathematics & Cybernetics.
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
