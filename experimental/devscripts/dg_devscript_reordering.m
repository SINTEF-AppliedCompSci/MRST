mrstModule add dg vem vemmech ad-props ad-core ad-blackoil blackoil-sequential gasinjection reorder matlab_bgl

%%

gravity reset off;

n = 10;
l = 1000;
G = computeGeometry(cartGrid([n,n], [l,l]*meter));
G = computeVEMGeometry(G);
G = computeCellDimensions(G);

rock = makeRock(G, 100*milli*darcy, 1);
fluid = initSimpleADIFluid('phases', 'WO'                   , ...
                           'rho'   , [1, 1]*kilogram/meter^3, ...
                           'mu'    , [0.5, 0.5]*centi*poise     , ...
                           'n'     , [1, 1]                 );

modelfi = TwoPhaseOilWaterModel(G, rock, fluid);
modelFV = getSequentialModelFromFI(modelfi);
modelDG = modelFV;

%%
jt = 0.6;
ot = 1e-3;
mt = 1e-3;
degree = 2;
disc   = DGDiscretization(modelDG.transportModel, ...
                         'degree', degree, ...
                         'basis' , 'legendre', ...
                         'useUnstructCubature', true,  ...
                         'jumpTolerance', jt, ...
                         'outTolerance', ot, ...
                         'meanTolerance', mt, ...
                         'plotLimiterProgress', false);
modelDG.transportModel = TransportOilWaterModelDG(G, rock, fluid, 'disc', disc, 'dsMaxAbs', 0.2);    
modelDG.pressureModel  = PressureOilWaterModelSemiDG(G, rock, fluid, 'disc', disc);    

%%

time = 1*year;
rate = 1*sum(poreVolume(G, rock))/time;
W = [];
W = addWell(W, G, rock, 1          , 'type', 'rate', 'val', rate    , 'comp_i', [1,0]);
W = addWell(W, G, rock, G.cells.num, 'type', 'bhp' , 'val', 50*barsa, 'comp_i', [1,0]);

dt    = 7*day;
dtvec = rampupTimesteps(time, dt, 0);
schedule = simpleSchedule(dtvec, 'W', W);

sW     = 0.0;
state0 = initResSol(G, 100*barsa, [sW,1-sW]);
state0 = assignDofFromState(modelDG.transportModel.disc, state0);

%%

[wsDG, statesDG, rep] = simulateScheduleAD(state0, modelDG, schedule);

%%
    
[modelDG.transportModel.extraStateOutput, modelDG.pressureModel.extraStateOutput] = deal(true);
modelDGreorder = modelDG;
modelDGreorder.pressureModel.extraStateOutput = true;

modelDGreorder.transportModel = ReorderingModelDG(modelDGreorder.transportModel, 'plotAfterCellSolve', false);

modelDGreorder.transportModel.chunkSize = 20;
modelDGreorder.transportModel.parent.extraStateOutput = true;


[wsDGReorder, statesDGReorder, repDGReorder] = simulateScheduleAD(state0, modelDGreorder, schedule);

%%

close all

figure
plotToolbar(G, statesDG); colormap(jet);
figure


plotToolbar(G, statesDGReorder); colormap(jet);

plotWellSols({wsDG, wsDGReorder});

%%

mrstModule add upr
n = 5;
G = pebiGrid(l/n, [l,l]);
G = computeVEMGeometry(G);
G = computeCellDimensions(G);

rock = makeRock(G, 100*milli*darcy, 1);

W = [];
xw = [0,l/2; l,l/2];
isInj = [true, false];
for wNo = 1:size(xw,1)
    d = sqrt(sum((xw(wNo,:) - G.cells.centroids).^2, 2));
    wc = find(d == min(d));
    wc = wc(1);
    if isInj(wNo)
        W = addWell(W, G, rock, wc, 'type', 'rate', 'val', rate, 'comp_i', [1,0]);
    else
        W = addWell(W, G, rock, wc, 'type', 'bhp', 'val', 50*barsa, 'comp_i', [1,0]);
    end
end

dt    = 30*day;
dtvec = rampupTimesteps(time, dt, 0);
schedule = simpleSchedule(dtvec, 'W', W);

%%
    
modelfi = TwoPhaseOilWaterModel(G, rock, fluid);
modelFV = getSequentialModelFromFI(modelfi);
modelDG = modelFV;

%%

[jt, ot, mt] = deal(Inf);
jt = 0.2;
ot = 1e-3;
degree = 0;
disc   = DGDiscretization(modelDG.transportModel, 'degree', degree, ...
                         'basis', 'legendre', 'useUnstructCubature', true,  'jumpTolerance', jt, ...
     'outTolerance', ot, 'meanTolerance', mt);
modelDG.transportModel = TransportOilWaterModelDG(G, rock, fluid, 'disc', disc);


sW     = 0.0;
state0 = initResSol(G, 100*barsa, [sW,1-sW]);
state0 = assignDofFromState(modelDG.transportModel.disc, state0);

%%

[modelDG.transportModel.extraStateOutput, modelDG.pressureModel.extraStateOutput] = deal(true);
modelDGreorder = modelDG;
modelDGreorder.pressureModel.extraStateOutput = true;

modelDGreorder.transportModel = ReorderingModelDG_ghost(modelDGreorder.transportModel, 'plotProgress', false);

modelDGreorder.transportModel.chunkSize = 1;
modelDGreorder.transportModel.parent.extraStateOutput = true;


[wsDGReorder, statesDGReorder, repDGReorder] = simulateScheduleAD(state0, modelDGreorder, schedule);

%%

[wsDG, statesDG, rep] = simulateScheduleAD(state0, modelDG, schedule);

%%

close all

figure
plotToolbar(G, statesDG); colormap(jet);
figure
plotToolbar(G, statesDGReorder); colormap(jet);
plotWellSols({wsDG, wsDGReorder});

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
