%% 1D Tutorial For a Oil-Water-Surfactant system
% The input data is read from a deck using Eclipse format
% (SURFACTANT1D.DATA). The surfactant property (see file surfact.inc) are taken
% from SPE paper 145036.
%
% Surfactant is added to water in order to decrease the surface tension so that,
% in particular, the residual oil is mobilized. See more detail about the
% modeling equations in ad-eor/docs
%
% In this example, water and surfactant are injected at the left-hand side and
% oil is produced at the right-hand side at a given pressure.
%
% In a first period, only water is injected. Then, for a second period,
% surfactant is added to the water.

%% We load the necessary modules
%

mrstModule add ad-core ad-blackoil ad-eor ad-props deckformat mrst-gui


%% We load the input data and setup the grid, rock and fluid structures
% 

current_dir = fileparts(mfilename('fullpath'));
if mrstIsLiveEditorDir(current_dir)
   % Running in "Live Editor" cell mode.  Fall back to expected location.

   current_dir = fullfile(mrstPath('ad-eor'), 'examples', 'surfactant');
end

fn = fullfile(current_dir, 'SURFACTANT1D.DATA');
gravity off

deck = readEclipseDeck(fn);
deck = convertDeckUnits(deck);

fluid = initDeckADIFluid(deck);
G = initEclipseGrid(deck);
G = computeGeometry(G);
rock  = initEclipseRock(deck);
rock  = compressRock(rock, G.cells.indexMap);


%% Set up the initial state
% Constant pressure, residual water saturation, no surfactant
%

nc = G.cells.num;
state0 = initResSol(G, 280*barsa, [ .2, .8]); % residual water saturation is 0.2
state0.cs    = zeros(G.cells.num, 1);
state0.csmax = state0.cs;

%% Set up the model
% 
% The model object contains the grid, the fluid and rock properties and the
% modeling equations. See simulatorWorkFlowExample.
%


model = OilWaterSurfactantModel(G, rock, fluid, ...
                                                  'inputdata', deck, ...
                                                  'extraStateOutput', true);                                      

%% Convert the deck schedule into a MRST schedule by parsing the wells
%

schedule = convertDeckScheduleToMRST(model, deck);

%% Run the schedule
%
% We use the function simulateScheduleAD to run the simulation
% Options such as maximum non-linear iterations and tolerance can be set in
% the system struct.

[wellSolsSurfactant, statesSurfactant] = simulateScheduleAD(state0, model, ...
                                                  schedule);

scheduleOW = schedule;
scheduleOW.control(2).W(1).cs = 0;
scheduleOW.control(2).W(2).cs = 0;
                                              
[wellSolsOW, statesOW] = simulateScheduleAD(state0, model, ...
                                                  scheduleOW);

figure()
plotToolbar(G, statesSurfactant, 'startplayback', true, 'plot1d', true, 'field', 's:1');

plotWellSols({wellSolsSurfactant, wellSolsOW}, ...
             cumsum(schedule.step.val), ...
             'DatasetNames', {'Surfactant', 'Water flooding'});

%% Copyright notice
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
