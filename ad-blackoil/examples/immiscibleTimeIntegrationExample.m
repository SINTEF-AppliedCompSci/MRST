%% Notice on Computational Cost
warning('ComputationalCost:High', ...
       ['Please be advised that this example often takes a long time ', ...
        'to run (e.g., more than 30 minutes of CPU time)']);
pause(10)

%%
mrstModule add ad-core ad-blackoil deckformat ad-props linearsolvers
mrstVerbose on
L = 1000*meter;
G = cartGrid(150, L);
G = computeGeometry(G);

x = G.cells.centroids/L;

poro = repmat(0.5, G.cells.num, 1);
%%
A = [0.20, 0.30, 0.05];
B = [0.45, 0.55, 0.1];
C = [0.70, 0.80, 0.05];
regs = {A, B, C};
nreg = numel(regs);
for i = 1:nreg
    R = regs{i};
    poro(x > R(1) & x < R(2)) = R(3);
end
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
implicit = packSimulationProblem(state0, model, schedule, 'immiscible_time', 'Name', 'FIM', 'Description', 'Fully-implicit');
%% Explicit solver
% This solver has a time-step restriction based on the CFL condition in
% each cell. The solver estimates the time-step before each solution.
model_explicit = setTimeDiscretization(model, 'explicit', 'verbose', 2);
explicit = packSimulationProblem(state0, model_explicit, schedule, 'immiscible_time',...
                                        'Name', 'EXPL', 'Description', 'Explicit');

fsb = model_explicit.FlowDiscretization.getFlowStateBuilder();
disp(fsb)


%% Adaptive implicit
% We can solve some cells implicitly and some cells explicitly based on the
% local CFL conditions. For many cases, this amounts to an explicit
% treatment far away from wells or driving forces. The values for
% estimated composition CFL and saturation CFL to trigger a switch to
% implicit status can be adjusted.
model_aim = setTimeDiscretization(model, 'adaptive-implicit', 'verbose', 2);
aim = packSimulationProblem(state0, model_aim, schedule, 'immiscible_time', ...
            'Name', 'AIM', 'Description', 'AIM');
%% Make an explicit solver with larger CFL limit
% Since the equation is linear, we can set the NonLinearSolver to use a
% single step. We bypass the convergence checks and can demonstrate the
% oscillations that result in taking a longer time-step than the stable
% limit.
model_explicit_largedt = setTimeDiscretization(model, 'explicit', 'verbose', true,...
    'saturationCFL', 5, ...
    'compositionCFL', inf); % Immiscible, saturation cfl is enough

model_explicit_largedt.stepFunctionIsLinear = true; 
explicit_largedt = packSimulationProblem(state0, model_explicit_largedt, schedule, ...
    'immiscible_time', 'Description', 'Explicit (CFL>1)', 'Name', 'EXPL_unstable');

%% Simulate the problems
problems = {implicit, explicit, aim, explicit_largedt};
simulatePackedProblem(problems, 'continueOnError', false);
%% Get output and plot the well results
% There are oscillations in the well curves. Increasing the CFL limit
% beyond unity will eventually lead to oscillations.
[ws, states, reports] = getMultiplePackedSimulatorOutputs(problems);
names = cellfun(@(x) x.Description, problems, 'UniformOutput', false);
%% Plot the results
% Note oscillations for CFL > 1.
model_cfl = model.setupStateFunctionGroupings();
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
    cfl = estimateSaturationCFL(model_cfl, state, dt);
    plotCellData(G, cfl);
    plotCellData(G, ones(G.cells.num, 1), 'Color', 'k');
    legend('Porosity', 'CFL', 'CFL stable limit');
    drawnow
end
%%
ns = numel(states);
figure(1);
for stepNo = 1:numel(schedule.step.val)
    clf; hold on
    for r = 1:nreg
        R = regs{r};
        x1 = L*R(1);
        x2 = L*R(2);
        patch([x1, x1, x2, x2], [0, 1, 1, 0], 'k', 'FaceAlpha', .1);
        text((x1 + x2)/2, 0.8, sprintf('\\phi = %1.2f', R(3)),'HorizontalAlignment', 'center', 'FontSize', 14)
    end
    h = nan(ns, 1);
    for j = 1:ns
        if j == 1
            st = '--';
            lw = 2;
        else
            st = '-';
            lw = 1;
        end
        h(j) = plotCellData(G, states{j}{stepNo}.s(:, 1), 'LineStyle', st, 'LineWidth', lw);
    end
    set(gca, 'FontSize', 13);
    legend(h, names);
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
%% Solve with sequential version of the solver
% In the previous section, we used the variable implicitness, but the
% solver did still simultaneously update pressure and saturations. An
% alternative is to use a sequential solver which first solves pressure,
% and then solves the transport equations either entirely or fully
% implicitly.
mrstModule add sequential
modelseq = getSequentialModelFromFI(model);

seq = packSimulationProblem(state0, modelseq, schedule, 'immiscible_time', 'Description', 'Sequential-Implicit', 'Name', 'SI');

modelseq_explicit = setTimeDiscretization(modelseq, 'explicit', 'verbose', true);
seq_explicit = packSimulationProblem(state0, modelseq_explicit, schedule, 'immiscible_time', 'Description', 'Sequential-Explicit', 'Name', 'SE');

modelseq_aim = setTimeDiscretization(modelseq, 'explicit', 'verbose', true);
seq_aim = packSimulationProblem(state0, modelseq_aim, schedule, 'immiscible_time', 'Description', 'Sequential-AIM', 'Name', 'SAIM');

allvariants = {implicit, explicit, aim, seq, seq_explicit, seq_aim};
simulatePackedProblem(allvariants, 'continueOnError', false);

[ws, states, reports, snames, T] = getMultiplePackedSimulatorOutputs(allvariants);
names = cellfun(@(x) x.Description, problems, 'UniformOutput', false);

%% Plot the displacement
% Note that as this example has no changes in mobility or significant
% compressibility, we can expect the different levels of explicitness to be
% the same between the sequential and fully coupled solvers.
markers = {'-.', '.', 'o'};
for i = 1:numel(T{1})
    clf; hold on
    for j = 1:numel(states)
        plot(states{j}{i}.s(:, 1), markers{mod(j, 3)+1});
    end
    legend(names)
    drawnow
end
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
