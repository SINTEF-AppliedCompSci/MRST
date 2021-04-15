gravity off
mrstModule add ntpfa ad-core ad-blackoil ad-props ...
               sequential mimetic incomp vem upr

dims = [15, 15];
pdims = [1, 1];

G = cartGrid(dims, pdims);
G = twister(G, 0.1);
G = computeGeometry(G);
fluid = initSimpleADIFluid();

rock = makeRock(G, 1, 1);
model = PressureOilWaterModelNTPFA(G, rock, fluid);
model.verbose = true;

[bc, src, W] = deal([]);
bc = pside(bc, G, 'xmin', 1, 'sat', [1, 0]);
bc = pside(bc, G, 'xmax', 0, 'sat', [1, 0]);
% bc = [];
% src = addSource(src, 1, 1, 'sat', [1, 0]);
% src = addSource(src, G.cells.num, -1, 'sat', [1, 0]);

%%
% NTPFA
state0 = initResSol(G, 0, [1, 0]);
dT = 1;

schedule = simpleSchedule(dT, 'W', W, 'bc', bc, 'src', src);
[~, states] = simulateScheduleAD(state0, model, schedule);
ntpfa = states{1};
% MIMETIC
S = computeMimeticIP(G, rock);
fluid2 = initSimpleFluid('mu', [1, 1], 'rho', [1, 1], 'n', [1, 1]);
mimetic = incompMimetic(state0, G, S, fluid2, 'bc', bc, 'src', src);

% TPFA
T = computeTrans(G, rock);
fluid2 = initSimpleFluid('mu', [1, 1], 'rho', [1, 1], 'n', [1, 1]);
tpfa = incompTPFA(state0, G, T, fluid2, 'bc', bc, 'src', src);
%%
figure(1), clf
subplot(1, 3, 1);
plotCellData(G, ntpfa.pressure)
title('N-TPFA')

subplot(1, 3, 2);
plotCellData(G, mimetic.pressure)
title('Mimetic')

subplot(1, 3, 3);
plotCellData(G, tpfa.pressure)
title('TPFA')


figure; hold on
x = G.cells.centroids(:, 1);
plot(x, ntpfa.pressure, '.')
plot(x, mimetic.pressure, '.')
plot(x, tpfa.pressure, '.')
legend('NTPFA', 'Mimetic', 'TPFA')

%%
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
