%% Diagnostics GUI using MRST simulation output
%
% This example shows how to use the diagnostics post processing GUI and
% visualise results from MRST simulation output.
% First we run the Egg model in MRST using the packSimulationProblem
% setup. Then we load the results into the GUI which calculates the
% diagnostics and displays them interactively.
%
% See also PostProcessDiagnosticsECLIPSE.m
%
%% Load relevant modules
mrstModule add ad-blackoil ad-core deckformat mrst-gui ad-props diagnostics example-suite

%% Get Eclipse deck for EGG
% Here we perform simulations on a single realization of the EGG model
deck = getDeckEGG('realization',1);
gravity reset on

%% Setup model, schedule and initial state
% We get all necessary parameters from the Eclipse deck
G = initEclipseGrid(deck);
[state0, model, schedule, nonlinear] = initEclipseProblemAD(deck, 'G', G);

%% Define/pack simulation problem
problem = packSimulationProblem(state0, model, schedule, 'egg_model_FlowDiagnostics');

%% Simulate problem
% This may take some time to simulate, however if the simulation has been
% completed already a new call to simulatePackedProblem will recognise
% that the simulation is complete and will do nothing.
% If the simulation is aborted, a new call to simulatePackedProblem will continue
% the simulation.
% To clear previous results and rerun a simulation use:
%       clearPackedSimulatorOutput(problem);
simulatePackedProblem(problem);

%% states/wellSols are accessed through handles
[ws, states, reports] = getPackedSimulatorOutput(problem);

%% Run PostProcessDiagnostics
% PostProcessDiagnosticsMRST takes as input a packed problem, packed using
% the packSimulationProblem function. The simulation must have been run
% already to get the flow field at each timestep.
% Diagnostics are calculated and displayed interactively.
PostProcessDiagnosticsMRST(problem);


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
