%% Example demonstrating the two-phase oil-water Egg model
% This example sets up and runs the Egg model using the two-phase AD
% solvers. 
%
% For details on the EggModel and the corresponding ensamble, see
% Jansen, J. D., et al. "The egg model–a geological ensemble for reservoir
% simulation." Geoscience Data Journal 1.2 (2014): 192-195.

mrstModule add ad-core ad-blackoil deckformat
% Realizations can be set to [] for base cae, or a number between 1 and 100
% for different permeabilities.
realization = [];
[G, rock, fluid, deck, state] = setupEGG('realization', realization);
model = selectModelFromDeck(G, rock, fluid, deck);
schedule = convertDeckScheduleToMRST(model, deck);

%% Run simulation
[wellSols, states] = simulateScheduleAD(state, model, schedule);
%% Plot the injector BHP and the well oil and water rates
% Since the injectors are rate controlled and the producers are pressure
% controlled, we can plot the quantities that vary. We also plot the water
% cut.

T = cumsum(schedule.step.val)/year;
W = schedule.control(1).W;

inj = find(vertcat(W.sign) > 0);
prod = find(vertcat(W.sign) < 0);
bhp = getWellOutput(wellSols, 'bhp', inj);

clf,
plot(T, bhp/barsa, 'linewidth', 2)
legend({W(inj).name})
xlabel('Time [years]')
ylabel('Injector bottom-hole pressure [bar]')

orat = getWellOutput(wellSols, 'qOs', prod);
wrat = getWellOutput(wellSols, 'qWs', prod);

figure;
plot(T, -orat*day, 'linewidth', 2)
legend({W(prod).name})
xlabel('Time [years]')
ylabel('Producer oil-rate [m^3/day]');

figure;
plot(T, -wrat*day, 'linewidth', 2)
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
plotWellSols(wellSols)

figure;
plotToolbar(G, states)
plotWell(G, W(inj), 'fontsize', 12, 'color', 'k')
plotWell(G, W(prod), 'fontsize', 12, 'color', 'r')

view(-50, 70);
axis tight off