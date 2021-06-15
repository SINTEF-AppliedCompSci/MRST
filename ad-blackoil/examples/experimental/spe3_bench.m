mrstModule add ad-blackoil ad-core mrst-gui ad-unittest
[schedule, model, state0] = getBenchmarkAD('spe3');
G    = model.G;
rock = model.rock;

%%
% Plot fluid model
fluidPlotPanelAD(model);

% Plot permeability field
figure;
plotToolbar(G, rock)
title('Rock')
view(20, 60);
axis tight

% Plot initial state
figure;
plotToolbar(G, state0)
title('Rock')
view(20, 60);
axis tight

% Add well to figure
plotWell(G, schedule.control(1).W)

%% Simulate the schedule itself
[ws, states, report] = simulateScheduleAD(state0, model, schedule);
dt = schedule.step.val;
%% Simulate the schedule with dynamic timestepping
schedule_compressed = compressSchedule(schedule);
solver = getNonLinearSolver(model);
[ws_timesteps, states_timesteps, report_timesteps] = simulateScheduleAD(state0, model, schedule_compressed, ...
                            'outputministeps', true, 'nonlinearsolver', solver);

[~, dt_timesteps] = convertReportToSchedule(report_timesteps, schedule_compressed);
%%
% Plot well curves
plotWellSols({ws, ws_timesteps}, {dt, dt_timesteps});

% Plot reservoir properties during simulation
figure;
plotToolbar(G, states)
title('Simulation results')
view(20, 60);
axis tight
plotWell(G, schedule.control(1).W)

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
