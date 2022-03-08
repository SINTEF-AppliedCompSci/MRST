%% Enhanced geothermal system (EGS)
% This example shows how to simulate a small geothermal system in an
% artificial fracture newtork. The grid is constructed by extruding a 2D
% PEBI grid with refinement around volumetric fractures vertically.

%% Add necessary MRST modules
mrstModule add ad-core ad-props ad-blackoil
mrstModule add geothermal compositional
mrstModule add test-suite
mrstModule add upr
mrstModule add mrst-gui
mrstVerbose on

%% Set up example
test = TestCase('small_egs_geothermal');

%% Plot setup
test.figure();
plotGrid(test.model.G, test.model.G.cells.tag, 'faceColor', [1,1,1]*0.8, 'edgeColor', 'none');
plotGrid(test.model.G, 'faceColor', 'none', 'edgeAlpha', 0.1);
test.setAxisProperties(gca);
camlight()
plotWell(test.model.G, test.schedule.control(1).W, 'color', 'k', 'fontSize', 30);
axis off

%% Simulate
problem = test.getPackedSimulationProblem();
clearPackedSimulatorOutput(problem, 'prompt', true);
simulatePackedProblem(problem);

%% Interactive plot of results
close all
[wellSols, states, reports] = getPackedSimulatorOutput(problem);
test.plot(states);
a = 0.6; cmap = hot.*a + (1-a);
colormap(cmap);
plotWellSols(wellSols, test.schedule.step.val)

%% Temperature at selected timesteps
T = getWellOutput(wellSols, 'T');
Tmin = min(min(T));
Tmax = max(max(T));
steps = [8, 15, 30];
for i = 1:numel(steps)
   test.figure();
   plotCellData(test.model.G, states{steps(i)}.T, 'edgeAlpha', 0.2);
   test.setAxisProperties(gca);
   plotWell(test.model.G, test.schedule.control(1).W, 'color', 'k', 'fontSize', 0.01);
   colormap(cmap), camlight();
   axis off;
   caxis([Tmin, Tmax]);
end

%% Plot EGS efficiency
p   = getWellOutput(wellSols, 'bhp');
T   = getWellOutput(wellSols, 'T');
q   = getWellOutput(wellSols, 'qWs');

[h, rho] = deal(zeros(size(p)));
for i = 1:2
    h(:, i)   = test.model.fluid.hW(p(:,i), T(:,i));
    rho(:, i) = test.model.fluid.rhoW(p(:,i), T(:,i));
end
qH  = abs(q.*rho.*h);

eff = (qH(:,2))./qH(:,1);
time = cumsum(test.schedule.step.val);

figure('Position', [0, 0, 800, 200])

hold on
kWatt = kilo*joule/second;
plot(time/year, eff, 'color', 'k', 'linew', 2);
% Indicate timesteps plotted above by circles
plot(time(steps)/year, eff(steps), 'ok', 'linew', 2);
axis([[time(5), time(end)]/year, min(eff(5:end))*0.95, max(eff(5:end))*1.05]);
set(gca, 'Box', true, 'FontSize', 13);
xlabel('Time (years)')

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