%% Example involving significant gravity effects in a CO2 injection scenario
mrstModule add compositional ad-core ad-props mrst-gui sequential
gravity reset on
%% Essential set up
BaseName = '2d_water_co2_gravity';
dims = [50, 1, 20];
T = 273.15 + 30;
p = 50*barsa;
time = 5*year;
dt = 10*day;
cnames = {'Water', 'CarbonDioxide'};
e = 1e-3;
z_0 = [1-e, e]; % Initial composition
z   = [e, 1-e]; % Injection composition
%% Make grid and fluid model
G = cartGrid(dims, [100, 10, 10]);
G = computeGeometry(G);
rock = makeRock(G, [0.1, 0.1, 0.001]*darcy, 0.3); % Kh = 100*Kv

f_eos = TableCompositionalMixture(cnames);
ncomp = f_eos.getNumberOfComponents();
eos = EquationOfStateModel(G, f_eos);
% Comment out the following to switch to K-values for VLE prediction
% K_wat = @(p, varargin) 0*p + 100;
% K_gas = @(p, varargin) 0*p + 1e-3;
% K_values = {K_wat, K_gas};
% eos.equilibriumConstantFunctions = K_values;
[~, ~, ~, ~, ~, rhoO] = standaloneFlash(p, T, [0, 1], eos);
[~, ~, ~, ~, ~, ~, rhoG] = standaloneFlash(p, T, [1, 0], eos);

fluid = initSimpleADIFluid('phases', 'OG', 'blackoil', false, 'rho', [rhoO, rhoG]);
fluid.krG = @(S) S.^2;
fluid.krO = @(S) S.^2;
%% Use overall composition formulation
model = GenericOverallCompositionModel(G, rock, fluid, f_eos, 'water', false);
%% Set up driving forces
irate = 0.5*sum(model.operators.pv)/time;
% Inject in lower left corner, produce in upper right corner
W = [];
W = verticalWell(W, G, rock, 1, 1, dims(3), 'type', 'rate', 'val', irate, 'components', z, 'comp_i', [0, 1]);
W = verticalWell(W, G, rock, dims(1), 1, 1, 'type', 'bhp', 'val', p, 'components', z_0, 'compi', [0, 1]);

rhoS = [rhoO, rhoG]; % Set surface densities in well at res conditions
W(1).rhoS = rhoS;
W(2).rhoS = rhoS;
%% Set up schedule (timesteps + forces)
timesteps = rampupTimesteps(time, dt);
schedule = simpleSchedule(timesteps, 'W', W);
state0 = initCompositionalState(G, p, T, [1, 0], z_0, eos);
%% Simulate the base case with fully-implicit phase-potential upwind
nls = NonLinearSolver('useRelaxation', true);
fim = packSimulationProblem(state0, model, schedule, BaseName, 'NonLinearSolver', nls, 'Name', 'Fully-implicit PPU');
simulatePackedProblem(fim, 'restartStep', 1, 'continueOnError', false);
%% Plot the gas plume
[~, states] = getPackedSimulatorOutput(fim);
%%
figure;
plotToolbar(G, states, 'field', 's:2');
plotGrid(G, vertcat(W.cells), 'FaceColor', 'none', 'EdgeColor', 'r');
view(0, 0);
%% Solve the regular phase-upwind
% The Brenier & Jaffre approach is used to correctly phase-upwind the
% quantities as the mobility changes.
pmodel = PressureModel(model);
tmodel = TransportModel(model);
smodel = SequentialPressureTransportModel(pmodel, tmodel);
% Use residual-based pressure convergence criteria
smodel.pressureModel.useIncTol   = false;
smodel.pressureModel.pressureTol = 1e-3;
seq = packSimulationProblem(state0, smodel, schedule, BaseName, 'NonLinearSolver', nls, 'Name', 'Sequential-PU');
simulatePackedProblem(seq, 'restartStep', 1, 'continueOnError', false);
%% Solve with hybrid upwind that treats gravity separately from pressure gradient ("viscous" flow)
smodel_hu = setHybridUpwindDiscretization(smodel);
seq_hu = packSimulationProblem(state0, smodel_hu, schedule, BaseName, 'NonLinearSolver', nls, 'Name', 'Sequential-HU');
simulatePackedProblem(seq_hu, 'restartStep', 1, 'continueOnError', false);
%%
problems = {fim, seq, seq_hu};
[~, states, reports, names] = getMultiplePackedSimulatorOutputs(problems);
n = numel(states);
%% Plot the transport iterations together with the solver output
stats = cellfun(@(x) getPressureTransportIterations(x), reports, 'UniformOutput', false);
t = value(cellfun(@(x) x.transport, stats', 'UniformOutput', false));
t = cumsum(t);
nstep = numel(schedule.step.val);
for i = 1:nstep
    figure(1); clf;
    for solverNo = 1:n
        subplot(2, n, solverNo)
        plotCellData(G, states{solverNo}{i}.s(:, 2), 'EdgeColor', 'none');
        title(names{solverNo});
        view(0, 0);
        axis tight equal
        daspect([1, 1, 0.2])
    end
    subplot(2, n, solverNo+1:n*2);
    plot(t(1:i, 2:end));
    legend(names(2:end));
    xlim([0, nstep]);
    ylim([0, 1.1*max(t(:))]);
    xlabel('Step');
    ylabel('Cumulative transport iterations')
    drawnow
end
%% Plot state functions graphs for transport
figure;
plotStateFunctionGroupings(model, 'Stop', 'ComponentTotalFlux');
title(names{1});
figure;
plotStateFunctionGroupings(smodel.transportModel, 'Stop', 'ComponentTotalFlux');
title(names{2});
figure;
plotStateFunctionGroupings(smodel_hu.transportModel, 'Stop', 'ComponentTotalFlux');
title(names{3});
%% Show the functions used
for i = 2:3
    tmp = problems{i}.SimulatorSetup.model.validateModel();
    disp(names{i});
    tmp.transportModel.parentModel.FlowDiscretization
end
%% Copyright Notice
%
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
