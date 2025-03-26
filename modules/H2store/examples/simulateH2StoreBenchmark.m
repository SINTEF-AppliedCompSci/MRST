%% Modified 3D Benchmark of Hydrogen Storage Cycling in an Aquifer
% This example adapts a benchmarking simulation originally designed to test gas conversion
% in gas storage systems, focusing on cyclic hydrogen (H₂) storage in an aquifer.
% original Benchmark:
% Hogeweg, S., Strobel, G., & Hagemann, B. (2022). Benchmark study for the simulation of underground
% hydrogen storage operations. Comput Geosci, 26, 1367–1378.
% https://doi.org/10.1007/s10596-022-10163-5
%
% For further details on this case, refer to:
% Ahmed, E., et al. (2024). Phase behavior and black-oil simulations of hydrogen storage in saline aquifers. 
% Advances in Water Resources, 191, 104772.
%
% Key modifications include:
% - The initial state is filled with water, incorporating hydrostatic boundary pressure 
%   due to the low compressibility of water.
% - The black-oil model from MRST is employed, where oil represents the water phase, and 
%   the water phase itself is disabled.
% - Hydrogen is injected into an initially brine-filled reservoir from an injection well 
%   located at the top of the reservoir.

% This model enables a detailed simulation of the interactions between hydrogen storage cycles, 
% the reservoir's behavior, and the impacts of various operational strategies.

% Additional characteristics:
% - Porosity and permeability are heterogeneously distributed, with an average porosity of 
%   15% and a mean horizontal permeability of 143 mD (vertical permeability approximately 
%   3 mD).
% - The Brooks–Corey model for capillary pressure is utilized, with parameters 
%   𝜆 = 2.0 and 𝑃𝑒 = 0.1 bar.
% - The simulation retains a conversion phase as a build-up phase, contrasting with the 
%   original setup where the reservoir was initially saturated with a gas mixture; 
%   our initial state is brine.

% The build-up phase consists of six cycles, each lasting 120 days of injection, followed by 
% a one-month idle period. After the build-up phase, we explore a cyclic injection/production 
% scenario that extends over 10 cycles.

clearvars;
%% Notice on Computational Cost
warning('ComputationalCost:High', ...
       ['Please be advised that this example takes a long time ', ...
        'to run']);
%% Necessary MRST modules for simulation
mrstModule add ad-core ad-blackoil ad-props deckformat mrst-gui

% Define the name of the simulation for output files
name = 'H2_STORAGE_RS_UHS_BENCHMARK';
%% Read the Eclipse deck file containing the simulation data
% Change input file by UHS_BENCHMARK_RS_SALT.DATA for SALT EFFECTS
baseDir = fileparts(mfilename('fullpath'));
dataFile = fullfile(baseDir,  'data', 'UHSModifiedBenchmark', 'UHS_BENCHMARK_RS.DATA');
deck = readEclipseDeck(dataFile);
% Note that the user can generate different PVT tables (see thermodynamics)

%% Prepare simulation parameters and initial state
[~, options, state0, model, schedule, ~] = ModelForUHSModifiedBenchmark(deck);

%% Add custom output functions to the model for additional diagnostics
model.OutputStateFunctions{end + 1} = 'CapillaryPressure';  % Output capillary pressure
model.OutputStateFunctions{end + 1} = 'SurfaceDensity';     % Output surface density
model.OutputStateFunctions{end + 1} = 'ShrinkageFactors';   % Output shrinkage factors
model.outputFluxes = false;                                 % Disable output of fluxes

%% Set up the linear and nonlinear solvers
lsolve = selectLinearSolverAD(model);                          % Select the linear solver for the model
nls = NonLinearSolver();                                       % Create a nonlinear solver object
nls.LinearSolver = lsolve;                                     % Assign the linear solver to the nonlinear solver

%% Pack the simulation problem with the initial state, model, and schedule
problem = packSimulationProblem(state0, model, schedule, name, 'NonLinearSolver', nls);

%% Run the simulation
tic;
simulatePackedProblem(problem,'restartStep',1);
elapsedTime = toc;
fprintf('Simulation completed in %.2f seconds. \n', elapsedTime);
%% gGet reservoir and well states
[ws,states] = getPackedSimulatorOutput(problem);
%% Plot well states
plotWellSols(ws);
%%
figure; plotToolbar(model.G, states); title(['MRST ', name]);
%% Compute total gas volume in the reservoir
nstep = numel(schedule.step.val);
dt = schedule.step.val;
time = cumsum(dt);
timeSteps = (1:nstep)';
[totGV_all, vH2InW_all, vH2InG_all] = computeGasVolumes(model, states, timeSteps);

%% Initialize variables for hydrogen volume, pressure, and rates
hydrogenLossPercentage = zeros(1, nstep);
totalGasVolume = zeros(1, nstep);
avPressure_res = zeros(1, nstep);
avSat_res = zeros(1, nstep);
bhp_res = zeros(1, nstep);
qGs_res = zeros(1, nstep);
qWs_res = zeros(1, nstep);

%% Compute cumulative hydrogen volume and averaged pressure over time
for step = 1:nstep
    totalGasVolume(step) = sum(totGV_all{step});
    dissolvedHydrogenVolume = sum(vH2InW_all{step});
    avPressure_res(step) = mean(states{step}.pressure);
    avSat_res(step) = mean(states{step}.s(:,2));
    
    if totalGasVolume(step) > 0
        hydrogenLossPercentage(step) = (dissolvedHydrogenVolume / totalGasVolume(step)) * 100;
    else
        hydrogenLossPercentage(step) = 0;
    end
    
    bhp_res(step) = ws{step}.bhp;
    qGs_res(step) = ws{step}.qGs;
    qWs_res(step) = ws{step}.qOs;
end

%% Plotting results over time (in months)
timeMonths = convertTo(time, 30.*day);

figure;
subplot(4,1,1);
plot(timeMonths, avSat_res, 'LineWidth', 1.5);
ylabel('Average Saturation');
title('Average Saturation Over Time');

subplot(4,1,2);
plot(timeMonths, avPressure_res, 'LineWidth', 1.5);
ylabel('Average Pressure (Pa)');
title('Average Pressure Over Time');

subplot(4,1,3);
plot(timeMonths, totalGasVolume, 'LineWidth', 1.5);
ylabel('Total Gas Volume');
title('Total Gas Volume Over Time');

subplot(4,1,4);
plot(timeMonths, cellfun(@sum, vH2InW_all), 'LineWidth', 1.5);
ylabel('Dissolved H2 Volume');
title('Dissolved H2 Volume Over Time');
xlabel('Time (months)');

figure;
subplot(4,1,1);
plot(timeMonths, bhp_res, 'LineWidth', 1.5);
ylabel('BHP (Pa)');
title('Bottom Hole Pressure Over Time');

subplot(4,1,2);
plot(timeMonths, qGs_res, 'LineWidth', 1.5);
ylabel('Gas Rate');
title('Gas Rate Over Time');

subplot(4,1,3);
plot(timeMonths, qWs_res, 'LineWidth', 1.5);
ylabel('Water Rate');
title('Water Rate Over Time');

subplot(4,1,4);
plot(timeMonths, hydrogenLossPercentage, 'LineWidth', 1.5);
ylabel('H2 Loss (%)');
title('Hydrogen Loss Percentage Over Time');
xlabel('Time (months)');
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


