%%
mrstModule add ad-props deckformat mrst-gui ad-core ad-blackoil
[G, rock, fluid, deck, state] = setupSPE1();
model = selectModelFromDeck(G, rock, fluid, deck);


schedule = convertDeckScheduleToMRST(model, deck);

schedule.step.val = schedule.step.val(1:5);
schedule.step.control = schedule.step.control(1:5);

model.extraStateOutput = true;
%% Run the entire schedule
[wellSols, states, report] = simulateScheduleAD(state, model, schedule);
%% Output surface fluxes from the stored properties
states = convertReservoirFluxesToSurface(model, states);

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
