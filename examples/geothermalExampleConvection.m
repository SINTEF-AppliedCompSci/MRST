%% Geothermal convection
% This example shows how you can simulate geothermal convection using
% sources and boundary conditions

%% Add modules
mrstModule add geothermal compositional
mrstModule add ad-core ad-props
mrstModule add upr
mrstModule add test-suite
mrstModule add mrst-gui

%% Setup 1: heating for 3000 years
setup1  = TestCase('convection_geothermal');
problem1 = setup1.getPackedSimulationProblem();
simulatePackedProblem(problem1, 'restartStep', 1);
[~, states1] = getPackedSimulatorOutput(problem1);

%% Set visaulization properties
Gviz = setup1.model.G;
Gviz.griddim = 2;
Gviz.cells.centroids = Gviz.cells.centroids(:,[1,3]);
Gviz.nodes.coords = Gviz.nodes.coords(:, [1,3]);
cmap = hot(100); a = 0.5; cmap = cmap.*a + (1-a);
cax  = [min(states1{1}.T), max(setup1.schedule.control(1).src.T)];

%% Visualize setup 1
fig1 = setup1.figure(); colormap(cmap); caxis(cax);
for i = 1:numel(states1)
    set(0, 'CurrentFigure', fig1), cla; hold on
    unstructuredContour(Gviz, states1{i}.T, 'fill', true);
    axis equal tight; set(gca, 'YDir', 'reverse', 'YLim', [-100, 4000]); box on
    v = faceFlux2cellVelocity(setup1.model.G, states1{i}.flux);
    v = v(:, [1,3]);
    quiver(Gviz.cells.centroids(:,1), Gviz.cells.centroids(:,2), v(:,1), v(:,2))
    drawnow(), pause(0.05);
end

%% Setup 2: heat for the first 500 years only
setup2   = TestCase('convection_geothermal',  'timeHeat', 500*year);
problem2 = setup2.getPackedSimulationProblem();
simulatePackedProblem(problem2, 'restartStep', 1);
[~, states2] = getPackedSimulatorOutput(problem2);

%% Visualize setup 2
fig2 = setup2.figure(); colormap(cmap); caxis(cax);
for i = 1:numel(states2)
    set(0, 'CurrentFigure', fig2), cla; hold on
    unstructuredContour(Gviz, states2{i}.T, 'fill', true);
    axis equal tight; set(gca, 'YDir', 'reverse', 'YLim', [-100, 4000]); box on
    v = faceFlux2cellVelocity(setup2.model.G, states2{i}.flux);
    v = v(:, [1,3]);
    quiver(Gviz.cells.centroids(:,1), Gviz.cells.centroids(:,2), v(:,1), v(:,2))
    drawnow(), pause(0.05);
end

%% Setup 3: heat for the first 500 years only, inject mass during heating
setup3   = TestCase('convection_geothermal', 'timeHeat', 500*year, 'massInjection', true);
problem3 = setup3.getPackedSimulationProblem();
simulatePackedProblem(problem3, 'restartStep', 1);
[~, states3] = getPackedSimulatorOutput(problem3);

%% Visualize setup 3
fig3 = setup3.figure(); colormap(cmap); caxis(cax);
for i = 1:numel(states3)
    set(0, 'CurrentFigure', fig3), cla; hold on
    unstructuredContour(Gviz, states3{i}.T, 'fill', true);
    axis equal tight; set(gca, 'YDir', 'reverse', 'YLim', [-100, 4000]); box on
    v = faceFlux2cellVelocity(setup1.model.G, states3{i}.flux);
    v = v(:, [1,3]);
    quiver(Gviz.cells.centroids(:,1), Gviz.cells.centroids(:,2), v(:,1), v(:,2))
    drawnow(), pause(0.05);
end

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
