%% Bucley-Leverett displacement
% In this example, we consider an incompressible Buckley-Leverett
% displacement in a horizontal 1D channel aligned with the x-axis, and
% compare higher-order discontinuous Galerkin methods with the standard
% finite-volume discretization.

mrstModule add dg ad-core ad-props ad-blackoil sequential
saveeps = @(a,b) disp(b);  % Dummy function

%% Set up problem
% Construct grid, compute geometry and cell dimensions, and set
% petrophysical properties                               
G    = computeGeometry(cartGrid([30,1]));
G    = computeCellDimensions(G);
rock = makeRock(G, 1, 1);

% We consider Bucley-Leverett-type displacement with quadratic relperms
fluid = initSimpleADIFluid('phases', 'WO' , 'n', [2,2], ... 
                           'mu', [1,1],  'rho', [1,1]);

% The base model is a generic black-oil model with oil and water
model  = GenericBlackOilModel(G, rock, fluid, 'gas', false);

% Initial state: filled with oil and unit volumetric flow rate
state0 = initResSol(G, 1, [0,1]);
state0.flux(1:G.cells.num+1) = 1;

% Boundary conditions
bc = fluxside([], G, 'left' ,  1, 'sat', [1,0]); % Inflow
bc = fluxside(bc, G, 'right', -1, 'sat', [1,0]); % Outflow

% Define schedule: unit time steps, simulate almost to breakthrough
schedule = simpleSchedule(ones(30,1), 'bc', bc);

%% Finite volume simulation
% Use standard discretization in MRST: single-point upwind (SPU)
tmodel           = TransportModel(model);
[~, stFV, repFV] = simulateScheduleAD(state0, tmodel, schedule);

%% Lowest-order dG simulation
% Simulate using dG(0). This should be equivalent to SPU.
% The model is built as a wrapper around the base model, so that all
% physics are consistently handled by this. The base model is stored as the
% property ''parentModel''
tmodelDG0 = TransportModelDG(model, 'degree', 0);
disp(tmodelDG0)

%% Simulate dG(0)
[~, stDG0, repDG0] = simulateScheduleAD(state0, tmodelDG0, schedule);

%% First-order dG simulation
% TransportModelDG also takes optional arguments to DGDiscretization, such
% as degree. DGDiscretization supports using different polynomial order in
% each coordinate direction. We are simulating a 1D problem on a 2D grid,
% so we use higher order in the x-direction only.
tmodelDG1 = TransportModelDG(model, 'degree', [1,0]);

%% Limiters
% Higher-order dG methods typically exhibit spurious oscillations near
% discontinuities. We countermand this with a limiter that effectively
% introduces enough numerical diffusion to dampen oscillations. We use two
% limiters by default: a TVB (Total Variation Bounded) limiter that adjusts
% the gradient whenever the jump across an interface is greater than 0, and
% a physical limiter that scales the solution so that the minimum and
% maximum in each cell is within physical limits. To see the limiters in
% action, we use ''afterStepFn'' to plot the unlimited and limited solution
% after each timestep
tmodelDG1.storeUnlimited = true; % Store the unlimited state in each step
fn = plotLimiter(tmodelDG1, 'plot1d', true, 'n', 500); % afterStepFn

%% Simulate dG(1)
[~, stDG1, repDG1] = ...
    simulateScheduleAD(state0, tmodelDG1, schedule, 'afterStepFn', fn);

%% Higher-order dG simulations
% In principle, the dg module supports arbitrary order, as long as you can
% provide a sufficiently accurate quadrature rule. These are stored in
% e.g., ''getELEMENTCubaturePointsAndWeights'', ELEMENT \in {Triangle,
% Square, Tetrahedron, Cube}. We simulate the problem with dG(2)
tmodelDG2 = TransportModelDG(model, 'degree', [2,0]);
tmodelDG2.storeUnlimited = true;
fn = plotLimiter(tmodelDG2, 'plot1d', true, 'n', 500);

%% Simulate dG(2)
[~, stDG2, repDG2] = ...
    simulateScheduleAD(state0, tmodelDG2, schedule, 'afterStepFn', fn);
                                                    
%% Compare results
% Finally, we plot the results, and compare them to the exact solution
close all
% Plotting options
coords = getPlotCoordinates(G, 'n', 500, 'plot1d', true);    % Plot coords
opt    = {'linewidth', 1.5, 'coords', coords, 'plot1d', true}; % Plot options
s      = buckleyLeverettProfile('nW', 2, 'nN', 2);       % Exact sol
tsteps = [5,10,20]; % Timesteps

for t = tsteps
    figure('Position', [100, 100, 800, 400])
    % Tweak FV solution so we can plot it as dG
    stFV_plot      = stDG0{t};
    stFV_plot.sdof = stFV{t}.s;
    hold on
    plotSaturationDG(tmodelDG0.discretization, stDG0{t}, opt{:});
    plotSaturationDG(tmodelDG1.discretization, stDG1{t}, opt{:});
    plotSaturationDG(tmodelDG2.discretization, stDG2{t}, opt{:});
    plotSaturationDG(tmodelDG0.discretization, stFV_plot, ...
        'plot1d', true, 'coords', coords, 'linewidth', 3, 'lineStyle', '--');
    plot(coords.points(:,1), s(coords.points(:,1), t), 'linewidth', 1.5, 'color', 'k')
    hold off
    legend({'dG(0)', 'dG(1)', 'dG(2)', 'SPU', 'exact'});
    box on
    ylim([-0.2, 1.2])
end

%% Plot dG(1) solution
close all
pos = [0,0,800,400];

figure('Position', pos)
t = 15;
hold on
    plotSaturationDG(tmodelDG1.discretization, stDG1{t}.ul, opt{:}, 'linewidth', 1.5, 'color', 'k');
    plotSaturationDG(tmodelDG1.discretization, stDG1{t}, opt{:}, 'linewidth', 4, 'color', 'k', 'linestyle', '--');
hold off
box on;
ylim([-0.2, 1.2])
ax = gca;
ax.FontSize = 14;
bx = [19.5, 22.5, -0.1, 0.25];
rectangle('Position', [bx([1,3]), bx([2,4]) - bx([1,3])], 'linestyle', '--');
drawnow(); pause(0.5), saveeps('buckley-leverett', 'limiter');

figure('Position', [0,0,350,400])
hold on
    plotSaturationDG(tmodelDG1.discretization, stDG1{t}.ul, opt{:}, 'linewidth', 1.5, 'color', 'k');
    plotSaturationDG(tmodelDG1.discretization, stDG1{t}, opt{:}, 'linewidth', 4, 'color', 'k', 'linestyle', '--');
hold off
box on;
axis(bx)
ax = gca;
ax.FontSize = 14;
[ax.XTick, ax.YTick] = deal([]);
legend({'Unlimited', 'Limited'}, 'fontSize', 15, 'Location', 'northeast')
drawnow(); pause(0.5), saveeps('buckley-leverett', 'limiter-zoom');

%% Compare FV and dG(0) - dG(2)

t = 15;
figure('Position', [100, 100, 800, 400])
% Tweak FV solution so we can plot it as dG
stFV_plot      = stDG0{t};
stFV_plot.sdof = stFV{t}.s;
hold on
plotSaturationDG(tmodelDG0.discretization, stDG0{t}, opt{:});
plotSaturationDG(tmodelDG1.discretization, stDG1{t}, opt{:});
plotSaturationDG(tmodelDG2.discretization, stDG2{t}, opt{:});
plotSaturationDG(tmodelDG0.discretization, stFV_plot, ...
    'plot1d', true, 'coords', coords, 'linewidth', 4, 'lineStyle', '--');
plot(coords.points(:,1), s(coords.points(:,1), t), 'linewidth', 1.5, 'color', 'k')
hold off
legend({'dG(0)', 'dG(1)', 'dG(2)', 'SPU', 'exact'}, 'fontsize', 14, 'location', 'southwest');
box on
ylim([-0.2, 1.2])
ax = gca;
ax.FontSize = 14;
bx = [16, 20, 0.35, 0.75];
rectangle('Position', [bx([1,3]), bx([2,4]) - bx([1,3])], 'linestyle', '--');
drawnow(); pause(0.5), saveeps('buckley-leverett', 'solutions')

figure('Position', [100, 100, 400, 400])
% Tweak FV solution so we can plot it as dG
stFV_plot      = stDG0{t};
stFV_plot.sdof = stFV{t}.s;
hold on
plotSaturationDG(tmodelDG0.discretization, stDG0{t}, opt{:});
plotSaturationDG(tmodelDG1.discretization, stDG1{t}, opt{:});
plotSaturationDG(tmodelDG2.discretization, stDG2{t}, opt{:});
plotSaturationDG(tmodelDG0.discretization, stFV_plot, ...
    'plot1d', true, 'coords', coords, 'linewidth', 4, 'lineStyle', '--');
plot(coords.points(:,1), s(coords.points(:,1), t), 'linewidth', 1.5, 'color', 'k')
hold off
box on
axis(bx)
ax = gca;
[ax.XTick, ax.YTick] = deal([]);
saveeps('buckley-leverett', 'solutions-zoom')

%%
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
