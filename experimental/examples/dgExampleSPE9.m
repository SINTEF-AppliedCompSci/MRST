mrstModule add dg vem vemmech ad-props ad-core ad-blackoil blackoil-sequential gasinjection mrst-gui

%%

[G, rock, fluid, deck, state0] = setupSPE9();

[ii, jj, kk] = gridLogicalIndices(G);
kkmax = 5;

keep = kk <= kkmax;
G = removeCells(G, ~keep);
G = computeGeometry(G);
G = computeVEMGeometry(G);
G = computeCellDimensions(G);
rock = compressRock(rock, keep);

% Determine the model automatically from the deck. It should be a
% three-phase black oil model with gas dissoluton.
model = selectModelFromDeck(G, rock, fluid, deck);

f = model.fluid;
pref = 100*barsa;
fluid = initSimpleADIFluid('phases', 'WO', ...
                           'n', [1,1], ...
                           'mu', [f.muW(pref), f.muG(pref)], ...
                           'rho', [f.rhoWS, f.rhoOS], ...
                           'c', [1e-10, 1e-5]/barsa);

modelFV = TwoPhaseOilWaterModel(model.G, model.rock, fluid);

% Convert the deck schedule into a MRST schedule by parsing the wells
schedule = convertDeckScheduleToMRST(modelFV, deck);

modelFV = getSequentialModelFromFI(modelFV);

W = schedule.control(1).W;

W = [];

time = 2*year;
rate = 0.25*sum(poreVolume(G, rock))/time;
bhp = 50*barsa;
W = verticalWell(W, G, rock, 1,1,[], 'type', 'rate', 'val', rate, 'comp_i', [1,0]);
W = verticalWell(W, G, rock, 24,25,[], 'type', 'bhp', 'val', bhp, 'comp_i', [1,0]);

dt = 30*day;
dtvec = rampupTimesteps(time, dt);
schedule = simpleSchedule(dtvec, 'W', W);
% [G, rock] = deal(modelFV.G, modelFV.rock);


modelFV.G = G;

degree = 0;
dNo = 1;
jt = 1;
ot = Inf;
mt = Inf;

modelDG = getSequentialModelFromFI(modelFV);
% modelDG.AutoDiffBackend = DiagonalAutoDiffBackend();
disc = DGDiscretization(modelDG.transportModel, 3, ...
        'degree', degree(dNo), ...
        'basis', 'legendre', ...
        'useUnstructCubature', true, ...
        'jumpTolerance', jt, ...
        'outTolerance', ot, ...
        'meanTolerance', mt);
modelDG.transportModel = TransportOilWaterModelDG(G, rock, fluid, 'disc', disc);    
% modelDG.transportModel.AutoDiffBackend = DiagonalAutoDiffBackend();
state0.s = state0.s(keep, [1,2]);
state0.pressure = state0.pressure(keep);

%%

state0 = assignDofFromState(modelDG.transportModel.disc, state0);

[wellsols, states, reports] =...
    simulateScheduleAD(state0, modelDG, schedule);

%%

[wellsolsFV, statesFV, reportsFV] =...
    simulateScheduleAD(state0, modelFV, schedule);

%% Copyright Notice
%
% <html>
% <p><font size="-1">
% Copyright 2009-2020 SINTEF Digital, Mathematics & Cybernetics.
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
