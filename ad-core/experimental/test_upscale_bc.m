mrstModule add ad-props ad-core ad-blackoil
G = cartGrid([10, 10, 10]);
G = computeGeometry(G);

W = [];
src = [];
bc = [];
bc = pside(bc, G, 'xmax', 100*barsa, 'sat', [1, 0]);
bc = fluxside(bc, G, 'xmin', 100, 'sat', [1, 0]);

src = addSource(src, 1:G.cells.num, 1, 'sat', [1, 0]);
schedule = simpleSchedule(1, 'bc', bc, 'W', W, 'src', src);

p = partitionUI(G, [2, 2, 2]);
CG = generateCoarseGrid(G, p);

rock = makeRock(G, 1, 1);
fluid = initSimpleADIFluid();

model = TwoPhaseOilWaterModel(G, rock, fluid);
model_c = upscaleModelTPFA(model, p);

schedule_c = upscaleSchedule(model_c, schedule, 'bcUpscaleMethod', 'linear');

% src = [];
% addSource(

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
