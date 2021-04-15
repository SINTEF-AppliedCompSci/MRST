mrstModule add dg ad-core ad-props ad-blackoil blackoil-sequential weno ...
    vem vemmech mrst-gui
mrstVerbose on

%%

[G, rock, fluid, deck, state0] = setupSPE1();
G = computeCellDimensions2(G);

model = selectModelFromDeck(G, rock, fluid, deck);

% Convert the deck schedule into a MRST schedule by parsing the wells
schedule = convertDeckScheduleToMRST(model, deck);

%%

modelSI = getSequentialModelFromFI(model);
% modelSI.transportModel.conserveWater = true;

%%

[wsFV, stFV, repFV] = simulateScheduleAD(state0, modelSI, schedule);

%%

transportModel = TransportBlackOilModelDG(G, rock, fluid, ...
                               'disgas', modelSI.transportModel.disgas, ...
                               'vapoil', modelSI.transportModel.vapoil, ...
                               'degree', 1);
modelDG        = SequentialPressureTransportModelDG(modelSI.pressureModel, transportModel);

[wsDG, stDG, repDG] = simulateScheduleAD(state0, modelDG, schedule);

%%

plotWellSols({wsDG, wsFV})

%%

sd = cellfun(@(s1, s2) compareStates(s1, s2), stFV, stDG, 'unif', false);
close all
plotToolbar(G, sd);

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
