%% Enhanced Geothermal System (EGS) Example
% This example demonstrates how to simulate a small enhanced geothermal system (EGS)
% in an artificial fracture network using MRST. The grid is constructed by extruding
% a 2D PEBI grid with vertical refinement around volumetric fractures. The workflow
% covers model setup, simulation, and post-processing of results, including efficiency analysis.
%
% Requirements:
%   - MRST core and the following modules: ad-core, ad-props, ad-blackoil, geothermal,
%     compositional, test-suite, upr, mrst-gui
%   - MATLAB R2021a or newer is recommended
%   - This script should be run from the MRST root or with MRST on the MATLAB path
%
% Output:
%   - Plots of the grid, well configuration, temperature evolution, and EGS efficiency
%   - Console output indicating simulation progress

%% Add necessary MRST modules
% Load core and required modules for geothermal simulation and visualization.
% The 'mrstModule add' commands ensure all dependencies are available.
mrstModule add ad-core ad-props ad-blackoil
mrstModule add geothermal compositional
mrstModule add test-suite
mrstModule add upr
mrstModule add mrst-gui
mrstVerbose on

%% Set up example case
% Create a test case for a small EGS geothermal model. The TestCase object
% sets up the grid, rock, fluid, wells, and schedule for the simulation.
test = TestCase('small_egs_geothermal');

%% Plot initial grid and well configuration
% Visualize the grid and tag cells, overlaying the well locations.
% The first plot shows cell tags (e.g., fracture, matrix), the second overlays
% the grid structure and wells for reference. This helps verify the model setup.
test.figure();
plotGrid(test.model.G, test.model.G.cells.tag, 'faceColor', [1,1,1]*0.8, 'edgeColor', 'none');
plotGrid(test.model.G, 'faceColor', 'none', 'edgeAlpha', 0.1);
test.setAxisProperties(gca);
camlight();
plotWell(test.model.G, test.schedule.control(1).W, 'color', 'k', 'fontSize', 30);
axis off;

drawnow; % Ensure plots are rendered before simulation

%% Simulate geothermal system
% Pack the simulation problem and run the simulation. This step assembles
% the model, schedule, and initial state, then runs the time-stepping solver.
problem = test.getPackedSimulationProblem();
clearPackedSimulatorOutput(problem, 'prompt', true); % Clear previous results if any
simulatePackedProblem(problem); % Run the simulation

disp('Simulation complete. Proceeding to post-processing...');

%% Interactive plot of results
% Retrieve and plot simulation results interactively. This section visualizes
% the evolution of temperature and other state variables over time.
close all;
[wellSols, states, reports] = getPackedSimulatorOutput(problem);
test.plot(states); % Interactive plot of temperature field
% Adjust colormap for better temperature contrast
alpha = 0.6; cmap = hot.*alpha + (1-alpha);
colormap(cmap);
plotWellSols(wellSols, test.schedule.step.val); % Plot well rates/pressures

drawnow;

%% Plot temperature at selected timesteps
% Visualize temperature distribution at specific timesteps to observe thermal
% breakthrough and cooling of the reservoir. The steps array selects which
% timesteps to plot.
T = getWellOutput(wellSols, 'T');
Tmin = min(min(T));
Tmax = max(max(T));
steps = [8, 15, 30]; % Example: early, mid, and late timesteps
for i = 1:numel(steps)
   test.figure();
   plotCellData(test.model.G, states{steps(i)}.T, 'edgeAlpha', 0.2);
   test.setAxisProperties(gca);
   plotWell(test.model.G, test.schedule.control(1).W, 'color', 'k', 'fontSize', 0.01);
   colormap(cmap), camlight();
   axis off;
   caxis([Tmin, Tmax]);
   title(sprintf('Temperature at timestep %d', steps(i)));
end

drawnow;

%% Plot EGS efficiency over time
% Calculate and plot the thermal efficiency of the EGS system.
% Efficiency is defined as the ratio of produced to injected enthalpy.
% This metric indicates how effectively the system extracts heat.
p   = getWellOutput(wellSols, 'bhp'); % Bottom-hole pressure for each well
T   = getWellOutput(wellSols, 'T');   % Well temperature
q   = getWellOutput(wellSols, 'qWs'); % Well flow rates

[h, rho] = deal(zeros(size(p)));
for i = 1:2
    h(:, i)   = test.model.fluid.hW(p(:,i), T(:,i));   % Enthalpy
    rho(:, i) = test.model.fluid.rhoW(p(:,i), T(:,i));  % Density
end
qH  = abs(q.*rho.*h); % Enthalpy flow rate (W)

eff = (qH(:,2))./qH(:,1); % Efficiency: produced/injected enthalpy
time = cumsum(test.schedule.step.val); % Simulation time (seconds)

figure('Position', [0, 0, 800, 200]);
hold on;
kWatt = kilo*joule/second;
plot(time/year, eff, 'color', 'k', 'linew', 2);
% Indicate timesteps plotted above by circles
plot(time(steps)/year, eff(steps), 'ok', 'linew', 2);
axis([[time(5), time(end)]/year, min(eff(5:end))*0.95, max(eff(5:end))*1.05]);
set(gca, 'Box', true, 'FontSize', 13);
xlabel('Time (years)');
ylabel('Thermal efficiency (produced/injected enthalpy)');
title('EGS System Efficiency Over Time');


disp('Simulation and post-processing complete. See figures for results.');

%% Copyright Notice
%
% <html>
% <p><font size="-1">
% Copyright 2009-2024 SINTEF Digital, Mathematics & Cybernetics.
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