% 2D Three-Phase Surfactant-Polymer Injection Case
%
% This example contains a simple 4000 m-by-200 m-by-125 m reservoir
% described on 20-by-1-by-5 uniform Cartesian grid. One injection well is
% located in the bottom two layers and one production well is located in
% the top two layers. Hydrostatic equilibration is used for initialization.
%
% The polymer injection schedule follows a typical polymer waterflooding
% strategy. The flooding process begins with primary waterflooding for 1260
% days, followed by a polymer slug injected over 1700 days, and then
% switching back to water injection. The injection well is under rate
% control with target rate 1000 m3/day and upper limit of 450 bar on the
% bottom-hole pressure (bhp), whereas the production well is under pressure
% control with target bottom-home pressure 260 bar.

clc
clear
close all

mrstModule add ad-core ad-blackoil ad-eor ad-props deckformat mrst-gui

%% Set up model and initial conditions
% The data required for the example
% The following are all the files needed for this tutorial
% Two files are the data for the simulation of surfactant polymer flooding.
current_dir = fullfile(mrstPath('ad-eor'), 'examples', 'surfactant-polymer');
fn = fullfile(current_dir, 'SURFACTANTPOLYMER2D.DATA');
gravity reset on;

deck = readEclipseDeck(fn);
% The deck is using metric system, MRST uses SI unit internally
deck = convertDeckUnits(deck);

% Construct physical model, initial state and dynamic well controls.
[state0, model, schedule] = ...
   initEclipseProblemAD(deck, 'UseLegacyModels', true);

model = model.validateModel();
figure
plotStateFunctionGroupings(model);

% Add initial surfactant & polymer concentration
state0.cp   = zeros([model.G.cells.num, 1]);
state0.cs   = zeros([model.G.cells.num, 1]);
% maximum surfactant & polymer concentration, used to handle the adsorption
state0.cpmax= zeros([model.G.cells.num, 1]);
state0.csmax= zeros([model.G.cells.num, 1]);

%% Select nonlinear and linear solvers

% Using physically normalized residuals for non-linear convergence
% calcuation.
model.useCNVConvergence = true;

% Setting up the non-linear solver.
nonlinearsolver = NonLinearSolver();
nonlinearsolver.useRelaxation = true;

%% Run the schedule with plotting function
% Once a system has been created it is trivial to run the schedule. Any
% options such as maximum non-linear iterations and tolerance can be set in
% the system struct.

% The AD-solvers allow for dyanmic plotting during the simulation process.
% We set up the following function to plot the evolution of the related
% variables (s:2 means oil saturation by default), the change of the well
% curves, and the a panel showing simulation progress. You can customize
% the function based on your own preference.
close all
fn = getPlotAfterStep(state0, model, schedule, ...
    'plotWell', true, 'plotReservoir', true, 'view', [20, 8], ...
                      'field', 's:2');
[wellSolsSP, statesSP, reportsSP] = ...
    simulateScheduleAD(state0, model, schedule, ...
                    'NonLinearSolver', nonlinearsolver, 'afterStepFn', fn);

% we use scheduleP to run the three phase black oil polymer flooding simulation.
scheduleP = schedule;
scheduleP.control(2).W(1).cs = 0;
fn1 = getPlotAfterStep(state0, model, scheduleP, ...
    'plotWell', true, 'plotReservoir', true, 'view', [20, 8], ...
                      'field', 's:2');
[wellSolsP, statesP, reportP] = simulateScheduleAD(state0, model, scheduleP, 'afterStepFn', fn1);

% we use scheduleS to run the three phase black oil surfactant flooding simulation.
scheduleS = schedule;
scheduleS.control(2).W(1).cp = 0;
fn2 = getPlotAfterStep(state0, model, scheduleS, ...
    'plotWell', true, 'plotReservoir', true, 'view', [20, 8], ...
                      'field', 's:2');
[wellSolsS, statesS, reportS] = simulateScheduleAD(state0, model, scheduleS, 'afterStepFn', fn2);

% we use scheduleW to run the three phase black oil water flooding simulation.
scheduleW = schedule;
scheduleW.control(2).W(1).cs = 0;
scheduleW.control(2).W(1).cp = 0;
fn3 = getPlotAfterStep(state0, model, scheduleW, ...
    'plotWell', true, 'plotReservoir', true, 'view', [20, 8], ...
                      'field', 's:2');
[wellSolsW, statesW, reportW] = simulateScheduleAD(state0, model, scheduleW, 'afterStepFn', fn3);
%% Plot cell oil saturation in different tsteps of differnt kinds of flooding

T = (80:23:268);

% Plot cell oil saturation in different tsteps of pure water flooding
sOmin = min( cellfun(@(x)min(x.s(:,2)), statesW) );
sOmax = max( cellfun(@(x)max(x.s(:,2)), statesW) );
figure
for i = 1 : length(T)
    subplot(3,3,i)
    plotCellData(model.G, statesW{T(i)}.s(:,2))
    plotWell(model.G, schedule.control(1).W, 'fontsize', 10)
    axis tight
    colormap(jet)
    view(3)
    caxis([sOmin, sOmax])
    title(['T = ', num2str(T(i))])
end
set(gcf, 'name', 'Oil saturation for water flooding')

% Plot cell oil saturation in different tsteps of polymer flooding
sOmin = min( cellfun(@(x)min(x.s(:,2)), statesP) );
sOmax = max( cellfun(@(x)max(x.s(:,2)), statesP) );
figure
for i = 1 : length(T)
    subplot(3,3,i)
    plotCellData(model.G, statesP{T(i)}.s(:,2))
    plotWell(model.G, schedule.control(1).W, 'fontsize', 10)
    axis tight
    colormap(jet)
    view(3)
    caxis([sOmin, sOmax])
    title(['T = ', num2str(T(i))])
end
set(gcf, 'name', 'Oil saturation for polymer flooding')

% Plot cell oil saturation in different tsteps of surfactant flooding
sOmin = min( cellfun(@(x)min(x.s(:,2)), statesS) );
sOmax = max( cellfun(@(x)max(x.s(:,2)), statesS) );
figure
for i = 1 : length(T)
    subplot(3,3,i)
    plotCellData(model.G, statesS{T(i)}.s(:,2))
    plotWell(model.G, schedule.control(1).W, 'fontsize', 10)
    axis tight
    colormap(jet)
    view(3)
    caxis([sOmin, sOmax])
    title(['T = ', num2str(T(i))])
end
set(gcf, 'name', 'Oil saturation for surfactant flooding')

% Plot cell oil saturation in different tsteps of surfactant-polymer flooding
sOmin = min( cellfun(@(x)min(x.s(:,2)), statesSP) );
sOmax = max( cellfun(@(x)max(x.s(:,2)), statesSP) );
figure
for i = 1 : length(T)
    subplot(3,3,i)
    plotCellData(model.G, statesSP{T(i)}.s(:,2))
    plotWell(model.G, schedule.control(1).W, 'fontsize', 10)
    axis tight
    colormap(jet)
    view(3)
    caxis([sOmin, sOmax])
    title(['T = ', num2str(T(i))])
end
set(gcf, 'name', 'Oil saturation for surfactant-polymer flooding')

%% Plot well solutions of water flooding and surfactant-polymer flooding
% The orange line denotes pure water flooding while the blue line denotes SP
% flooing

plotWellSols({wellSolsSP, wellSolsW},cumsum(schedule.step.val))

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
