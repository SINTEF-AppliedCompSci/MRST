%% Polymer solver with boundary conditions and sources
%
% This example is made just to illustrate how one can setup a polymer
% simulation with boundary conditions and/or source.
% 

% Required modules
mrstModule add deckformat ad-core ad-blackoil ad-eor ad-props mrst-gui


%% Setup case

% Grid, rock and fluid
deck  = readEclipseDeck('POLYMER.DATA');
deck  = convertDeckUnits(deck);
G     = initEclipseGrid(deck);
G     = computeGeometry(G);
rock  = initEclipseRock(deck);
rock  = compressRock(rock, G.cells.indexMap);
fluid = initDeckADIFluid(deck);

% Gravity
gravity on

% Initial state
state0      = initResSol(G, 100*barsa, [ .1, .9]);
state0.cp    = zeros(G.cells.num,1);
state0.cpmax = zeros(G.cells.num,1);

% Create model
model = OilWaterPolymerModel(G, rock, fluid);

% Setup some schedule
dt = 25*day;
nt = 40;
clear schedule
timesteps = repmat(dt, nt, 1);

%% Pressure (Dirichlet) Boundary Condition

% Create Dirichlet boundary condition
bc = pside([], G, 'xmin', 500*barsa, 'sat', [1 0]);
bc = pside(bc, G, 'xmax', 100*barsa, 'sat', [0 0]);
bc.cp = 4.*ones(size(bc.sat,1), 1);

schedule = simpleSchedule(timesteps, 'bc', bc);

% Simulate
[~, states] = simulateScheduleAD(state0, model, schedule);

% Plot results in GUI
figure;
plotToolbar(G, states,'field', 's:1','lockCaxis',true);
view([-10, 14]);
axis tight;
colorbar; caxis([0 1]);

%% Flux (Neumann) Boundary Condition

% Create Neumann boundary condition
bc = fluxside([], G, 'xmin',  0.005, 'sat', [1 0]);
bc = fluxside(bc, G, 'xmax', -0.005, 'sat', [0 0]);
bc.cp = 4.*ones(size(bc.sat,1), 1);
schedule = simpleSchedule(timesteps, 'bc', bc);

% Simulate
[~, states] = simulateScheduleAD(state0, model, schedule);

% Plot results in GUI
figure;
plotToolbar(G, states,'field', 's:1','lockCaxis',true);
view([-10, 14]);
axis tight;
colorbar; caxis([0 1]);


%% Source

% Create source
ijk = gridLogicalIndices(G);
srcCells = find(ijk{1}==5  & ijk{2}==5  & ijk{3}==2);
snkCells = find(ijk{1}==26 & ijk{2}==26 & ijk{3}==2);
srcVals  = 0.001.*ones(numel(srcCells),1);
src = addSource( [], srcCells,  srcVals, 'sat', [1 0]);
src = addSource(src, snkCells, -srcVals, 'sat', [0 0]);
src.cp = 4.*ones(size(src.sat,1), 1);
schedule = simpleSchedule(timesteps, 'src', src);

% Simulate
[~, states] = simulateScheduleAD(state0, model, schedule);

% Plot results in GUI
figure;
plotToolbar(G, states,'field', 's:1','lockCaxis',true);
view([-10, 14]);
axis tight;
colorbar; caxis([0 1]);

%% Copyright notice

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
