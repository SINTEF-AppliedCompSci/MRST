mrstModule add ad-core ad-blackoil deckformat ad-props linearsolvers
mrstVerbose on
L = 1000*meter;
G = cartGrid(150, L);
G = computeGeometry(G);

x = G.cells.centroids/L;

poro = repmat(0.5, G.cells.num, 1);
%%
poro(x > 0.2 & x < 0.3) = 0.05;
poro(x > 0.7 & x < 0.8) = 0.05;
poro(x > 0.45 & x < 0.55) = 0.1;

rock = makeRock(G, 1*darcy, poro);

close all
plot(poro)
%%
time = 500*day;

% Default oil-water fluid with unit values
fluid = initSimpleADIFluid('phases', 'WO', 'n', [1 1]);

% Set up model and initial state.
model = GenericBlackOilModel(G, rock, fluid, 'gas', false);
state0 = initResSol(G, 50*barsa, [0, 1]);

% Set up drive mechanism: constant rate at x=0, constant pressure at x=L
pv = poreVolume(G, rock);
injRate = sum(pv)/time;
bc = fluxside([], G, 'xmin', injRate, 'sat', [1, 0]);
bc = pside(bc, G, 'xmax', 0*barsa, 'sat', [0, 1]);

n  = 150;
dT = time/n;
timesteps = rampupTimesteps(time, dT, 8);
schedule = simpleSchedule(timesteps, 'bc', bc);

%% Fully-implicit
% The default discretization in MRST is fully-implicit. Consequently, we
% can use the model as-is.
implicit = packSimulationProblem(state0, model, schedule, 'immiscible_time', 'Name', 'Fully-Implicit');
%% Explicit solver
% This solver has a time-step restriction based on the CFL condition in
% each cell. The solver estimates the time-step before each solution.
model_explicit = setTimeDiscretization(model, 'explicit', 'verbose', true);
explicit = packSimulationProblem(state0, model_explicit, schedule, 'immiscible_time', 'Name', 'Explicit');

%% Adaptive implicit
% We can solve some cells implicitly and some cells explicitly based on the
% local CFL conditions. For many cases, this amounts to an explicit
% treatment far away from wells or driving forces. The values for
% estimated composition CFL and saturation CFL to trigger a switch to
% implicit status can be adjusted.
model_aim = setTimeDiscretization(model, 'adaptive-implicit', 'verbose', true);
aim = packSimulationProblem(state0, model_aim, schedule, 'immiscible_time', 'Name', 'Adaptive-Implicit (AIM)');
%% Make an explicit solver with larger CFL limit
% Since the equation is linear, we can set the NonLinearSolver to use a
% single step. We bypass the convergence checks and can demonstrate the
% oscillations that result in taking a longer time-step than the stable
% limit.
model_explicit_largedt = model.validateModel();
% Get the flux discretization
flux = model_explicit_largedt.FluxDiscretization;
% Use explicit form of flow state
fb = ExplicitFlowStateBuilder();
fb.saturationCFL = 5;
fb.compositionCFL = inf; % Immiscible, saturation cfl is enough
flux = flux.setFlowStateBuilder(fb);
model_explicit_largedt.FluxDiscretization = flux;
model_explicit_largedt.stepFunctionIsLinear = true; 
explicit_largedt = packSimulationProblem(state0, model_explicit_largedt, schedule, 'immiscible_time', 'Name', 'Explicit (CFL target 5)');

%% Simulate the problems
problems = {implicit, explicit, aim, explicit_largedt};
simulatePackedProblem(problems, 'continueOnError', false);
%% Get output and plot the well results
% There are oscillations in the well curves. Increasing the CFL limit
% beyond unity will eventually lead to oscillations.
[ws, states, reports, names, T] = getMultiplePackedSimulatorOutputs(problems);
%% Plot the results
% Note oscillations for CFL > 1.
figure(1);
for stepNo = 1:numel(schedule.step.val)
    clf; 
    subplot(2, 1, 1);
    hold on
    for j = 1:numel(states)
        plotCellData(G, states{j}{stepNo}.s(:, 1))
    end
    legend(names);
 
    subplot(2, 1, 2); hold on;
    plotCellData(G, rock.poro);
    dt = schedule.step.val(stepNo);
    state = states{1}{stepNo};
    cfl = estimateSaturationCFL(model, state, dt);
    plotCellData(G, cfl);
    plotCellData(G, ones(G.cells.num, 1), 'Color', 'k');
    legend('Porosity', 'CFL', 'CFL stable limit');
    drawnow
end
%% Plot the fraction of cells which are implicit
imp_frac = cellfun(@(x) sum(x.implicit)/numel(x.implicit), states{3});
figure;
plot(100*imp_frac);
title('Number of implicit cells in AIM')
ylabel('Percentage of implicit cells')
xlabel('Step index')
%% Plot the time-steps taken
% The implicit solvers use the control steps while the explicit solver must
% solve many substeps.
figure; hold on
markers = {'o', '.', 'x', '.'};
for stepNo = 1:numel(reports)
    dt = getReportMinisteps(reports{stepNo});
    x = (1:numel(dt))/numel(dt);
    plot(x, dt/day, markers{stepNo})
end
ylabel('Timestep length [day]')
xlabel('Fraction of total simulation time');
legend(names)

% <html>
% <p><font size="-1">
% Copyright 2009-2019 SINTEF Digital, Mathematics & Cybernetics.
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
