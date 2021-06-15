%% Tutorial for a simple 3D case of a oil-water-surfactant system
% This example contains a $31\times31\times3$ fine grid containing two injectors
% in opposite corners and one producer in the middle of the domain. All wells
% are completed in the top layers of cells.
%
% The schedule being used contains first a period of injection with surfactant,
% followed by a water flooding phase without surfactant. Finally, the water rate
% is reduced for the final time steps.
%
% The data is read from deck (SURFACTANT3D.DATA)
%

mrstModule add ad-core ad-blackoil ad-eor ad-props deckformat mrst-gui

%% We load the input data and setup the grid, rock and fluid structures
%

current_dir = fileparts(mfilename('fullpath'));
if mrstIsLiveEditorDir(current_dir)
   % Running in "Live Editor" cell mode.  Fall back to expected location.

   current_dir = fullfile(mrstPath('ad-eor'), 'examples', 'surfactant');
end

fn = fullfile(current_dir, 'SURFACTANT3D.DATA');
gravity on

deck = readEclipseDeck(fn);
deck = convertDeckUnits(deck);

fluid = initDeckADIFluid(deck);
G = initEclipseGrid(deck);
G = computeGeometry(G);
rock  = initEclipseRock(deck);
rock  = compressRock(rock, G.cells.indexMap);



%% Set up the initial state
% We want a layer of oil on top of the reservoir and water on the bottom.  To do
% this, we alter the initial state based on the logical height of each
% cell. Initially, the pressure is constant and the surfactant concentration is
% zero.

ijk = gridLogicalIndices(G);

state0 = initResSol(G, 300*barsa, [ .9, .1]);
state0.s(ijk{3} == 1, 2) = .9;
state0.s(ijk{3} == 2, 2) = .8;

% Enforce s_w + s_o = 1;
state0.s(:,1) = 1 - state0.s(:,2);

state0.cs    = zeros(G.cells.num, 1);
state0.csmax = state0.cs;

clf
plotCellData(G, state0.s(:,2));
plotGrid(G, 'facec', 'none')
title('Oil concentration')
axis tight off
view(70, 30);
colorbar;


%% Set up the model
% The model object contains the grid, the fluid and rock properties and the
% modeling equations. See simulatorWorkFlowExample.
%
model = OilWaterSurfactantModel(G, rock, fluid, 'inputdata', deck);



%% Setup the schedule
% The wells and wells controls obtained from the input deck are parsed and a
% MRST schedule is set up.

schedule = convertDeckScheduleToMRST(model, deck);

%% Run the schedule
%
% We use the function simulateScheduleAD to run the simulation
% Options such as maximum non-linear iterations and tolerance can be set in
% the system struct.
%
% Here, we also send a resulthandler to save the output data, see ResultHandler.
%

resulthandler = ResultHandler('dataDirectory', pwd, 'dataFolder', 'cache', 'cleardir', true);
[wellSolsSurfactant, statesSurfactant] = simulateScheduleAD(state0, model, ...
                                                  schedule, 'OutputHandler', ...
                                                  resulthandler);
plotToolbar(G,statesSurfactant,'field','s:1');
W = schedule.control(1).W;
view(70,30), plotWell(G, W), axis tight off

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
