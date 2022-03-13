mrstModule add dg ad-core ad-props ad-blackoil blackoil-sequential weno vem vemmech
mrstVerbose on

%%

[GF, rockF, fluid, deck, state0F] = setupSPE1();
modelF = selectModelFromDeck(GF, rockF, fluid, deck);
[ii, jj, kk] = gridLogicalIndices(GF);
keep = ii == 1 & kk == 1 & jj <= 3;
G = removeCells(GF, ~keep);
rock = compressRock(rockF, keep);
scheduleF = convertDeckScheduleToMRST(modelF, deck);

%%

G = computeCellDimensions2(G);
model = ThreePhaseBlackOilModel(G, rock, fluid, 'disgas', true);
W = [];
f = sum(poreVolume(G, rock))./sum(poreVolume(GF, rockF))*1;
schedule = scheduleF;
W = addWell(W, G, rock, 1, 'type', schedule.control.W(1).type, ...
                           'val', schedule.control.W(1).val*f, ...
                           'refDepth', schedule.control.W(1).refDepth, ...
                           'comp_i', schedule.control.W(1).compi);
W = addWell(W, G, rock, G.cells.num, 'type', schedule.control.W(2).type, ...
                           'val', schedule.control.W(2).val*f, ...
                           'refDepth', schedule.control.W(1).refDepth, ...
                           'comp_i', schedule.control.W(2).compi);

schedule.control(1).W = W;
                       
state0.pressure = state0F.pressure(keep);
state0.s = state0F.s(keep,:);
state0.rs = state0F.rs(keep,:);
                       
%%

modelSI = getSequentialModelFromFI(model);
% modelSI.transportModel.conserveWater = true;x
%%

ix = 1:80;
subschedule = schedule;
subschedule.step.val     = subschedule.step.val(ix);
subschedule.step.control = subschedule.step.control(ix);

%%

[wsFV, stFV, repFV] = simulateScheduleAD(state0, modelSI, subschedule);

%%

transportModel = TransportBlackOilModelDG(G, rock, fluid, ...
                               'disgas', modelSI.transportModel.disgas, ...
                               'vapoil', modelSI.transportModel.vapoil, ...
                               'degree', 0);
modelDG        = SequentialPressureTransportModelDG(modelSI.pressureModel, transportModel);

[wsDG, stDG, repDG] = simulateScheduleAD(state0, modelDG, subschedule);

%%

plotWellSols({wsDG, wsFV})

%%

close all
sd = cellfun(@(st1, st2) compareStates(st1, st2), stFV, stDG);
plotToolbar(G, sd, 'plot1d', true);

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
