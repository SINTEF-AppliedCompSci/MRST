%% Compare sequential solver to fully implicit, applied to the SPE1 problem
% This example simulates the first SPE benchmark using both sequential and
% fully implicit solvers. The problem is a gas injection with significant
% density and mobility changes between phases, so the assumption of a
% constant total velocity for the sequential scheme is violated, and there
% are differences in production profiles between the two schemes. We also
% test a third approach, wherein the sequential scheme revisits the
% pressure solution if the mobilities change significantly, in order to
% convergence to the fully implicit solution.

mrstModule add ad-core ad-blackoil ad-props sequential example-suite

%% Set up the initial simulation model
% We use the existing setupSPE1 routine to handle the setup of all
% parameters, which are then converted into a fully implicit model and a
% schedule.
[G, rock, fluid, deck, state] = setupSPE1();

model = selectModelFromDeck(G, rock, fluid, deck);
schedule = convertDeckScheduleToMRST(model, deck);

%% Run the entire fully implicit schedule
% We simulate the schedule with a fully implicit scheme, i.e. where we
% solve for both saturations and pressure simultanously.
[wsFIMP, statesFIMP, repFIMP] = simulateScheduleAD(state, model, schedule);

%% Set up and simulate the schedule using a sequential implicit scheme
% We convert the fully implicit model into a sequential model, a special
% model which contains submodels for both the pressure and transport. In
% this scheme, we first solve the pressure equation implicitly to obtain
% total volumetric fluxes at reservoir conditions, as well as total well
% rates. The pressure and fluxes are subsequently used as input for the
% transport scheme, which advects the saturations in the total velocity
% field using a fractional flow formulation.

seqModel = getSequentialModelFromFI(model);
[wsSeq, statesSeq, repSeq] = simulateScheduleAD(state, seqModel, schedule);

%% Simulate the schedule with the outer loop option enabled
% In the sequential implicit scheme, the transport equations are derived by
% assuming a fixed total velocity. For certain problems, this assumption is
% not reasonable, and the total velocity may be strongly coupled to the
% changes in saturation during a timestep. In this case, the problem is a
% gas injection scenario where the injected gas has a much higher mobility
% and much lower density, leading to significantly different results
% between the fully implicit and the sequential implicit schemes.
%
% One way to improve the results of the sequential implicit scheme is to
% employ an outer loop when needed. This amounts to revisiting the pressure
% equation after the transport has converged and checking if the pressure
% equation is still converged after the saturations have changed. If the
% pressure equation has a large residual after the transport, we then
% resolve the pressure with the new estimates for saturations.

% Make a copy of the model
seqModelOuter = seqModel;
% Disable the "stepFunctionIsLinear" option, which tells the
% NonLinearSolver that it is not sufficient to do a single pressure,
% transport step to obtain convergence. We set ds tolerance to 0.001
% (which is also the default) which is roughly equivalent to the maximum
% saturation error convergence criterion used in the fully implicit solver.
seqModelOuter.stepFunctionIsLinear = false;
seqModelOuter.incTolSaturation = 1.0000e-03;


[wsOuter, statesOuter, repOuter] = simulateScheduleAD(state, seqModelOuter, schedule);

%% Plot the results
% We plot the results for the three different temporal discretization
% schemes. Since the water is close to immobile and the injector and
% producers are operated on gas and oil rates respectively, we plot the gas
% production and the bottom hole pressure for each well.
%
% We clearly see that there are substantial differences in the well curves
% due the changes in total velocity. The sequential implicit scheme with
% the outer loop enabled is very close to the fully implicit solution and
% will get closer if the tolerances are tightened.
T = cumsum(schedule.step.val);

wellSols = {wsFIMP, wsSeq, wsOuter};
states = {statesFIMP, statesSeq, statesOuter};
names = {'Fully-implicit', 'Sequential-implicit', 'Sequential-outer'};
markers = {'-', '--', '.'};

close all;
for i = 1:numel(wellSols)
    figure(1); hold on
    qo = -getWellOutput(wellSols{i}, 'qGs', 2)*day;
    plot(T/day, qo, markers{i}, 'linewidth', 2)
    xlabel('Time [days]');
    ylabel('Gas production [m^3/day]');
    
    figure(2); hold on
    bhp = getWellOutput(wellSols{i}, 'bhp', 2)/barsa;
    plot(T/day, bhp, markers{i}, 'linewidth', 2)
    xlabel('Time [days]');
    ylabel('Producer bottom hole pressure [bar]');
    
    figure(3); hold on
    bhp = getWellOutput(wellSols{i}, 'bhp', 1)/barsa;
    plot(T/day, bhp, markers{i}, 'linewidth', 2)
    xlabel('Time [days]');
    ylabel('Injector bottom hole pressure [bar]');
end

% Add legends to the plots
for i = 1:3
    figure(i);
    legend(names, 'location', 'northwest')
end

% Plot the different saturations and pressures
figure;
for i = 1:numel(names)
    s = states{i};
    
    subplot(2, numel(names), i)
    plotCellData(model.G, s{end}.s(:, 3))
    caxis([0, 1])
    axis tight
    title(names{i})
    xlabel('Gas saturation')
    
    subplot(2, numel(names), i + numel(names))
    plotCellData(model.G, s{end}.pressure)
    axis tight
    xlabel('Pressure')
end
%% Set up interactive plotting
% We finish the example by launching interactive viewers for the well
% curves, as well as the different reservoir quantities.
mrstModule add mrst-gui
wsol = {wsSeq, wsFIMP, wsOuter};
wnames = {'Sequential', 'FIMP', 'Outer loop'};
states = {statesSeq, statesFIMP, statesOuter};

for i = 1:numel(states)
    figure(i), clf
    plotToolbar(G, states{i});
    title(wnames{i});
    axis tight
    view(20, 60);
    plotWell(G, schedule.control(1).W)
end

plotWellSols(wsol, T, 'datasetnames', wnames)


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
