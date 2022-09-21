%% Example Title
% This example demonstrate a 1D flooding example implemented based on
% boundary condition and generic surfactant-polymer model. In this case,
% only water and polymer are involved and polymer is injecting from one
% end (left) to the other end (right).
mrstModule add ad-core ad-blackoil ad-eor ad-props deckformat mrst-gui

%% Run fine-grid simulation
fn = '1D_FLOODING.DATA';
[state0, model, schedule_deck, nlsolver] = initEclipseProblemAD(fn);
model.fluid.mixPar = 0.85;

dt = schedule_deck.step.val;
bc = fluxside([], model.G, 'xmin', schedule_deck.control.W(1).val, 'sat', 1);
bc = pside(bc, model.G, 'xmax', schedule_deck.control.W(2).val, 'sat', 1);

% the cp should only apply to the left boundary
bc.cp = [schedule_deck.control.W(:).cp];
schedule = simpleSchedule(dt, 'bc', bc);

[~, states] = simulateScheduleAD(state0, model, schedule, ...
                                 'NonLinearSolver', nlsolver);

%% visualizing some results
plt = @(x) plot(model.G.cells.centroids(:,1), x, 'LineWidth', 2);

figure('Position', [50, 500, 900, 300])

subplot(1,3,1)
plt(states{end}.s(:,1)); title('water saturation')

subplot(1,3,2)
plt(states{end}.cp); title('polymer concentration (kg/m^3)')

subplot(1,3,3)
plt(convertTo(states{end}.pressure, barsa)); title('pressure (barsa)')

% interactively visualizing the running result
figure
plotToolbar(model.G, states, 'plot1d', true, 'field', 'cp');

%% Copyright Notice
%
% <html>
% <p><font size="-1">
% Copyright 2009-2022 SINTEF Digital, Mathematics & Cybernetics.
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
