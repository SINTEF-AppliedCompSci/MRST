%% Example Title
% This example demonstrate a 1D flooding example implemented based on
% boundary condition and generic surfactant-polymer model. In this case,
% only water and polymer are involved and polymer is injecting from one
% end (left) to the other end (right).
mrstModule add ad-core ad-blackoil ad-eor ad-props deckformat

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
[~, states] = simulateScheduleAD(state0, model, schedule, 'NonLinearSolver', nlsolver);

%% visualizing some results
f=figure(1);
f.Position=[50 500 900 300];
subplot(1,3,1), hold all
plot(model.G.cells.centroids(:,1),states{end}.s(:,1),'LineWidth',2); title('water saturation')
subplot(1,3,2), hold all
plot(model.G.cells.centroids(:,1),states{end}.cp(:,1),'LineWidth',2); title ('polymer concentration (m^3/kg)')
subplot(1,3,3), hold all
plot(model.G.cells.centroids(:,1),states{end}.pressure(:,1)/barsa,'LineWidth',2); title('pressure (barsa)')

% interactively visualizing the running result
figure(2)
plotToolbar(model.G, states, 'plot1d', true, 'field', 'cp');