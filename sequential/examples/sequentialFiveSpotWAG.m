%% Set up and run a five-spot problem with water-alternating-gas (WAG) drive
% One approach to hydrocarbon recovery is to inject gas or water. In the
% water-alternating-gas approach, wells change between injection of gas and
% water to improve the sweep efficiency.
%
% In this example, we set up a five-spot injection pattern where the
% injector controls alternate between water and gas. The primary purpose of
% the example is not to demonstrate an ideal or realistic
% WAG-configuration, but rather how to set up a problem with changing well
% controls and to demonstrate how the sequential solvers perform on
% problems with large changes in mobility and well controls.
%
% We begin by loading the required modules

mrstModule add ad-core ad-blackoil ad-props ...
               sequential spe10 mrst-gui

%% Set up grid and rock structure
% We define a 50x50x1 grid, spanning a 1 km by 1km domain. The porosity is
% assigned via a Gaussian field and a synthethic permeability is computed
% using a standard correlation.

cartDims = [50, 50, 1];
physDims = [1000, 1000, 1]*meter;
G = cartGrid(cartDims, physDims);
G = computeGeometry(G);
gravity reset off

rng(0);
poro = gaussianField(cartDims, [0.05, 0.3], 3, 8);
poro = poro(:);
perm = poro.^3 .* (5e-5)^2 ./ (0.81 * 72 * (1 - poro).^2);
% Build rock object
rock = makeRock(G, perm, poro);
figure(1); clf
plotCellData(G, rock.poro)
axis equal tight
colorbar
title('Porosity')

%% Set up wells and schedule
% The injection scenario corresponds to injecting one pore volume of gas
% and water over 30 years at standard conditions. We first place wells in
% the corners of the domain set to fixed injection rates as well as a
% bottom hole pressure producer in the middle of the domain.
pv = poreVolume(G, rock);
T = 30*year;
irate = sum(pv)/(T*4);
% Function handle for easily creating multiple injectors
makeInj = @(W, name, I, J, compi) verticalWell(W, G, rock, I, J, [],...
    'Name', name, 'radius', 5*inch, 'sign', 1, 'Type', 'rate',...
    'Val', irate, 'comp_i', compi);
W = [];
W = makeInj(W, 'I1', 1,           1,           []);
W = makeInj(W, 'I3', cartDims(1), cartDims(2), []);
W = makeInj(W, 'I4', 1,           cartDims(2), []);
W = makeInj(W, 'I2', cartDims(1), 1,           []);

I = ceil(cartDims(1)/2);
J = ceil(cartDims(2)/2);
% Producer
W = verticalWell(W, G, rock, I, J, [], 'Name', 'P1', 'radius', 5*inch, ...
    'Type', 'bhp', 'Val', 100*barsa, 'comp_i', [1, 1, 1]/3, 'Sign', -1);
% Create two copies of the wells: The first copy is set to water injection
% and the second copy to gas injection.
[W_water, W_gas] = deal(W);
for i = 1:numel(W)
    if W(i).sign < 0
        % Skip producer
        continue
    end
    W_water(i).compi = [1, 0, 0];
    W_gas(i).compi   = [0, 0, 1];
end
% Simulate 90 day time steps, with a gradual increase of steps in the
% beginning as the mobility changes rapidly when wells begin injecting for
% the first time.
dT_target = 90*day;
dt = rampupTimesteps(T, dT_target, 10);

% Set up a schedule with two different controls. In the first control, the
% injectors are set to the copy of the wells we made earlier where water
% was injected. In the second control, gas is injected instead.
schedule = struct();
schedule.control = [struct('W', W_water);... % Water control 1
                    struct('W', W_gas)]; % Gas control 2
% Set timesteps
schedule.step.val = dt;
% Alternate between the gas and water controls every 90 days
schedule.step.control = (mod(cumsum(dt), 2*dT_target) >= dT_target) + 1;

% Plot the changes in schedule control
figure(1), clf
ctrl = repmat(schedule.step.control', 2, 1);
x = repmat(cumsum(dt/year)', 2, 1);
y = repmat([0; 1], 1, size(ctrl, 2));
surf(x, y, ctrl)
colormap(jet)
view(0, 90)
axis equal tight
set(gca, 'YTick', []);
xlabel('Time [year]')
title('Control changes over time: Red for gas injection, blue for water')

%% Set up fluid and simulation model
% We set up a three-phase fluid with quadratic relative permeability
% curves, significant viscosity contrast between the phases and
% compressibility for the oil and gas phases.
fluid = initSimpleADIFluid('phases',    'WOG', ...
                           'rho',       [1000, 700, 250], ...
                           'n',         [2, 2, 2], ...
                           'c',         [0, 1e-4, 1e-3]/barsa, ...
                           'mu',        [1, 4, 0.25]*centi*poise ...
                           );
% Set up three-phase, immiscible model with fully implicit discretization
model = ThreePhaseBlackOilModel(G, rock, fluid, 'disgas', false, 'vapoil', false);
% Create sequential model from fully implicit model
seqModel = getSequentialModelFromFI(model);
% Set up initial reservoir at 100 bar pressure and completely oil filled.
state = initResSol(G, 100*barsa, [0, 1, 0]);

%% Simulate the schedule using the sequential solver
% We simulate the full schedule using simulateScheduleAD
[wsSeq, statesSeq] = simulateScheduleAD(state, seqModel, schedule);

%% Simulate the schedule using a fully implicit solver
% For comparison purposes, we also solve the fully implicit case
[wsFIMP, statesFIMP] = simulateScheduleAD(state, model, schedule);

%% Plot producer rates
% We plot the producer water, gas and oil rates. We can see the effect of
% the alternating controls directly on the production curves and we use a
% stairstep plot since the underlying data is not continuous. 
df = get(0, 'DefaultFigurePosition');
figure('position', df.*[1, 1, 2, 1])
names = {'qWs', 'qOs', 'qGs'};
hold on
t = cumsum(schedule.step.val)/day;
colors = lines(3);

% q = getWellOutput(wsSeq, 'qWs'
for i = 1:numel(names)
%     subplot(1, 3, 1)
    qSeq = -getWellOutput(wsSeq, names{i}, 5);
    qFIMP = -getWellOutput(wsFIMP, names{i}, 5);
    stairs(t, qSeq*day, '-o', 'color', colors(i, :));
    stairs(t, qFIMP*day, '.', 'color', colors(i, :), 'linewidth', 2);
end
legend('Water, sequential', 'Water FIMP', 'Oil, sequential', ...
       'Oil FIMP', 'Gas, sequential', 'Gas FIMP');
xlabel('Time [days]')
ylabel('Production at standard conditions [m^3/day]')
axis tight

%% Plot the injector bottom hole pressure
% The injector bottom hole pressure has a bigger discrepancy between the
% two schemes, since the injectors use the cell total mobility to compute
% the relation between completion flux and bottom hole pressure. The
% deviation between the schemes is a direct result of the lagging mobility
% values in the pressure equation.
clf; hold on
colors = lines(4);

for i = 1:4
    bhpSeq = getWellOutput(wsSeq, 'bhp', i)/barsa;
    bhpFIMP = getWellOutput(wsFIMP, 'bhp', i)/barsa;

    stairs(t, bhpSeq, '-o', 'color', colors(i, :));
    stairs(t, bhpFIMP, '.', 'color', colors(i, :), 'linewidth', 2);
end
xlabel('Time [days]')
ylabel('Injector bottom hole pressure [bar]')
axis tight

%% Launch interactive plotting
% Finally we launch interactive viewers for the well curves and reservoir
% quantities.
plotWellSols({wsSeq, wsFIMP}, cumsum(schedule.step.val), ...
    'datasetnames', {'Sequential', 'FIMP'}, 'linestyles', {'-', '-.'})

figure;
plotToolbar(G, statesSeq);
axis tight
plotWell(G, W);
title('Sequential implicit')

figure;
plotToolbar(G, statesFIMP);
axis tight
plotWell(G, W);
title('Fully-implicit')

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

