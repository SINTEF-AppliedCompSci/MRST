%% Example demonstrating the two-phase oil-water Egg model
% This example sets up and runs the Egg model using the two-phase AD
% solvers. 
%
% For details on the EggModel and the corresponding ensamble, see
% Jansen, J. D., et al. "The egg model–a geological ensemble for reservoir
% simulation." Geoscience Data Journal 1.2 (2014): 192-195.

mrstModule add ad-core ad-blackoil deckformat example-suite
% Realizations can be set to 0 for base cae, or a number between 1 and 100
% for different permeabilities.
realization = 0;
[G, rock, fluid, deck] = setupEGG('realization', realization);
[state, model, schedule, nonlinear] = initEclipseProblemAD(deck, 'G', G, 'TimestepStrategy', 'none');
%% Run simulation
[wellSols, states, report] = simulateScheduleAD(state, model, schedule, 'NonLinearSolver', nonlinear);
%% Plot the injector BHP and the well oil and water rates
% Since the injectors are rate controlled and the producers are pressure
% controlled, we can plot the quantities that vary. We also plot the water
% cut.

T = convertTo(cumsum(schedule.step.val), year);
W = schedule.control(1).W;

inj = find(vertcat(W.sign) > 0);
prod = find(vertcat(W.sign) < 0);
bhp = getWellOutput(wellSols, 'bhp', inj);

clf,
plot(T, convertTo(bhp, barsa), 'linewidth', 2)
legend({W(inj).name})
xlabel('Time [years]')
ylabel('Injector bottom-hole pressure [bar]')

orat = getWellOutput(wellSols, 'qOs', prod);
wrat = getWellOutput(wellSols, 'qWs', prod);

figure;
plot(T, convertTo(-orat, meter^3/day), 'linewidth', 2)
legend({W(prod).name})
xlabel('Time [years]')
ylabel('Producer oil-rate [m^3/day]');

figure;
plot(T, convertTo(-wrat, meter^3/day), 'linewidth', 2)
legend({W(prod).name})
xlabel('Time [years]')
ylabel('Producer water-rate [m^3/day]');

figure;
plot(T, abs(wrat)./abs(wrat + orat), 'linewidth', 2)
legend({W(prod).name})
xlabel('Time [years]')
ylabel('Water cut');
%% Plot the pressure and water saturation through the simulation
df = get(0, 'DefaultFigurePosition');
h = figure('Position', df.*[1, 1, 2.25, 1]);
for i = 1:numel(states)
    % Plot the pressure
    timestr = formatTimeRange(sum(schedule.step.val(1:i)));
    figure(h); clf
    subplot(1, 2, 1)
    plotCellData(G, states{i}.pressure/barsa, 'EdgeColor', 'none');
    plotWell(G, W(inj), 'fontsize', 12, 'color', 'k')
    plotWell(G, W(prod), 'fontsize', 12, 'color', 'r')
    title(['Reservoir pressure after ', timestr]);
    view(-50, 70);
    axis tight off
    colorbar
    
    % Plot water saturation
    subplot(1, 2, 2)
    plotCellData(G, states{i}.s(:, 1), 'EdgeColor', 'none');
    plotWell(G, W(inj), 'fontsize', 12, 'color', 'k')
    plotWell(G, W(prod), 'fontsize', 12, 'color', 'r')
    title(['Water saturation after ', timestr]);
    view(-50, 70);
    axis tight off
    colorbar
    caxis([0, 1])
    drawnow
end
%% Launch interactive plotting
mrstModule add mrst-gui
plotWellSols(wellSols, cumsum(schedule.step.val))

figure;
plotToolbar(G, states)
plotWell(G, W(inj), 'fontsize', 12, 'color', 'k')
plotWell(G, W(prod), 'fontsize', 12, 'color', 'r')

view(-50, 70);
axis tight off

%% Copyright Notice
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
