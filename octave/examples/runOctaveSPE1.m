%% Demonstrate black-oil solver for SPE1 problem
% This example demonstrates how the AD-OO simulators can be run in Octave, with 
% a few caveats. Please note that the performance of MRST in Octave is far from
% optimal and the execution time for even simple problem can be an order of
% magnitude longer than in a recent version of Matlab. 

% Load modules
mrstModule add ad-props deckformat mrst-gui ad-core ad-blackoil
%% Set up problem
% We use the Octave module which loads a few compatibility fixes required. Note 
% that this code was written for Octave 4.2.1 and the fixes may not be 
% applicable for older or newer versions.
if exist('OCTAVE_VERSION', 'builtin') ~= 0
    mrstModule add octave
    run set_reasonable_octave_defaults;
end

[G, rock, fluid, deck, state] = setupSPE1();
% Due to a difference in behavior, the deck reader will read in commented out
% time-steps. We overcome this by loading them directly from file.
deck.SCHEDULE.step.val = getTimestepsSPE1_octave();
% Uniform control
deck.SCHEDULE.step.control = ones(size(deck.SCHEDULE.step.val));
%% Instansiate model
model = selectModelFromDeck(G, rock, fluid, deck);

% Convert the deck schedule into a MRST schedule by parsing the wells
schedule = convertDeckScheduleToMRST(model, deck);

% Plot setup
figure;
plotCellData(G, convertTo(rock.perm(:,1), milli*darcy), 'EdgeColor', 'k');
plotWell(G, schedule.control.W, 'radius',.5); 
title('Permeability (mD)')
axis tight, view(35, 40), colorbar('SouthOutside');
drawnow
%% Run simulation
[wellSols, states, report] = simulateScheduleAD(state, model, schedule);

%% Plot gas saturation
% The 3D plotting in Octave works with most MRST routines. However, the
% graphical user interfaces in MRST does not presently work with MRST.
figure;
for i = 1:numel(states)
  plotCellData(G, states{i}.s(:, 3), 'EdgeColor', 'k');
  plotWell(G, schedule.control.W, 'radius',.5); 
  title(['Gas saturation at ', formatTimeRange(report.ReservoirTime(i))]);
  axis tight, view(35, 40), colorbar('SouthOutside');
  drawnow
end

%%
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
