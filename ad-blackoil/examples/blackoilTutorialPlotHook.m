%% Example demonstrating in-situ plotting capabilities in MRST-AD
% The AD-solvers allow for dynamic plotting during the simulation process.
% This example demonstrates this capability using "plotWellSols" and a
% panel showing simulation progress.
%
% We first set up a simulation model of SPE1 in the standard manner.
mrstModule add ad-core ad-blackoil ad-props mrst-gui deckformat example-suite

[G, rock, fluid, deck, state0] = setupSPE1();
model = selectModelFromDeck(G, rock, fluid, deck);
schedule = convertDeckScheduleToMRST(model, deck);

%% Set up plotting function
% We set up a function handle that takes in the current simulator
% variables, which will be run after each succesful control step. This
% function can also pause the simulation or change the maximum timestep as
% a proof of concept.
close all
fn = getPlotAfterStep(state0, model, schedule, 'plotWell', true,...
                                               'plotReservoir', false);
disp(fn)

%% Run the simulation with plotting function
linsolve = selectLinearSolverAD(model);
[wellSols, states, report] = simulateScheduleAD(state0, model, schedule, ...
    'Verbose', true, 'afterStepFn', fn, 'linearsolver', linsolve);

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
