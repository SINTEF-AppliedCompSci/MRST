%% Illustrative Case: Efficiency of Structural Trap for Hydrogen Storage
% 
% This example simulates the storage of hydrogen in a saline aquifer and 
% investigates the influence of physical properties such as H2 solubility, 
% salinity, and molecular diffusion on various factors, including the 
% hydrogen plume, reservoir pressure, caprock tightness, and the recoverability 
% and loss of hydrogen. The model features a 2D dome-shaped aquifer with 
% integrated caprock and bedrock layers, extending horizontally for 50 m 
% and vertically for 50 m.  It represents the first case discussed 
% in the paper: "Phase behavior and black-oil simulations of hydrogen 
% storage in saline aquifers" (Ahmed, E., et al., 2024).
%
% In this simulation, various scenarios are explored regarding the 
% behavior of hydrogen in saline aquifers. The reader is encouraged 
% to examine different input files for specific case studies:
% - **H2STORAGE_RS.DATA**: Models only dissolved gas in brine, with 
%   evaporated water interactions disabled.
% - **H2STORAGE_RS_SALT.DAT**: Represents hydrogen interacting with brine 
%   at a molality of 4 mol/kg.
%
% The upper boundary of the aquifer is determined by the function 
% F(x) = σ + r sin(πx), where σ is set to 25 and r to 5, creating a 
% trapping structure for efficient H2 storage. The containment of injected 
% hydrogen gas relies on sufficiently high entry (capillary) pressure and 
% the low permeability of the caprock and bedrock formations.
%
% Hydrostatic pressure boundary conditions are applied to the lateral 
% boundaries, calculated from a reference pressure p_r at depth z = 0 of 
% 40 bar. No-flux conditions are enforced at the top and bottom boundaries. 
% The initial state of the brine occupies the domain in hydrostatic 
% equilibrium matching the boundary conditions. Physical and fluid parameters 
% are derived from Oko-roafor et al. (2023) and summarized in a related table. 
% The simulation also tests the injection of H2 into both pure water and brine 
% with salinity of 5 mol kg−1.
%--------------------------------------------------------------------------

clearvars; 
mrstModule add ad-core ad-blackoil ad-props deckformat mrst-gui upr

%% Define the case name and read the Eclipse deck file
name = 'H2_STORAGE_RS';
%% Use H2STORAGE_RS_SALT.DATA for brine
baseDir = fileparts(mfilename('fullpath')); % Get directory of the script
dataFile = fullfile(baseDir,  'data', 'Simple2DAquifer', 'H2STORAGE_RS.DATA');
deck = readEclipseDeck(dataFile);

%% Notice on Computational Cost
warning('ComputationalCost:High', ...
       ['Please be advised that this is a 10 cycle injection example which often takes some time ', ...
        'to run: reduce cycles for example']);

%% Set up the simulation parameters and model components
[~, ~, state0, model, schedule, ~] = ModelForSimple2DAquifer(deck);

%% Plot Grid with Wells, Permeability, and Porosity
figure;
subplot(1, 2, 1); 
plotCellData(model.G, model.rock.poro);
hold on;
plotGrid(model.G, schedule.control(1).W.cells, 'LineStyle', 'none', 'FaceColor', 'red');
title('Porosity');
axis off tight; 
subplot(1, 2, 2); 
plotCellData(model.G, log10(model.rock.perm(:, 1)));
hold on;
plotGrid(model.G, schedule.control(1).W.cells, 'LineStyle', 'none', 'FaceColor', 'red');
title('Permeability (log10)');
axis off tight; 
sgtitle('2D Dome-Shaped Aquifer'); 

%% Add output functions to the model for various properties
model.OutputStateFunctions{end + 1} = 'CapillaryPressure';
model.OutputStateFunctions{end + 1} = 'SurfaceDensity';
model.OutputStateFunctions{end + 1} = 'ShrinkageFactors';
model.outputFluxes = false;

%% Initialize the nonlinear solver and select the linear solver
nls = NonLinearSolver(); 
lsolve = selectLinearSolverAD(model); 
nls.LinearSolver = lsolve;

%% Pack the simulation problem with the defined components
problem = packSimulationProblem(state0, model, schedule, name, 'NonLinearSolver', nls);

%% Execute the simulation of the packed problem
simulatePackedProblem(problem, 'restartStep',1);

%% Get packed reservoir and well states
[ws, states] = getPackedSimulatorOutput(problem);
%% Plot states
figure;
plotToolbar(model.G, states);

%% Plot well output
figure;
plotWellSols(ws);

%% Get last steps of each injection, production, and idle period
LastSteps = getCyclesLastSteps(schedule);

%% Compute gas volumes for different phases
[totGV_charge, vH2InW_charge, vH2InG_charge] = computeGasVolumes(model, states, LastSteps.charge);
[totGV_disc, vH2InW_disc, vH2InG_disc] = computeGasVolumes(model, states, LastSteps.discharge);
[totGV_cush, vH2InW_cush, vH2InG_cush] = computeGasVolumes(model, states, LastSteps.cushion);

%% Add plots for H2 volume in Water, Gas, and Total Gas Volumes
figure; 
subplot(2, 2, 1);
plotCellData(model.G, states{LastSteps.cushion(1)}.s(:, 2));
title('s_G');
colorbar;
subplot(2, 2, 2);
plotCellData(model.G, vH2InG_cush{1});
title('V_{H2} in Gas');
colorbar;
subplot(2, 2, 3);
plotCellData(model.G, vH2InW_cush{1});
title('V_{H2} in Water');
colorbar;
subplot(2, 2, 4); 
plotCellData(model.G, totGV_cush{1});
title('Total V_{H2}');
colorbar;
sgtitle('V_{H2} after build-up phase'); 
%% Compute Hydrogen Loss Due to Dissolution in Percentage
month = 30*day;
Time = cumsum(schedule.step.val)/month;
hydrogenLossPercentage = zeros(1, length(LastSteps.discharge));
% Loop through each step of the discharge cycle
for step = 1:length(LastSteps.discharge)
    totalGasVolume = sum(totGV_disc{step});  % Total gas volume in the reservoir
    dissolvedHydrogenVolume = sum(vH2InW_disc{step});  % Volume of hydrogen dissolved in water
    if totalGasVolume > 0  % Check to avoid division by zero
        hydrogenLossPercentage(step) = (dissolvedHydrogenVolume / totalGasVolume) * 100;  % Percentage calculation
    else
        hydrogenLossPercentage(step) = 0;  % Set to 0 if no gas volume is present
    end
end
% Plot the hydrogen loss percentage over the discharge steps
figure;  % Create a new figure
plot(Time(LastSteps.discharge), hydrogenLossPercentage, 'LineWidth', 2, 'Color', 'b');
xlabel('Discharge Steps');
ylabel('Hydrogen Loss Percentage (%)');
title('Hydrogen Loss Due to Dissolution');
grid on;

%% Plot xH2 Over Time
rhoOS = model.fluid.rhoOS;  % Density of oil-saturated phase
rhoGS = model.fluid.rhoGS;  % Density of gas-saturated phase
StepsOfCycles = sort([LastSteps.cushion, LastSteps.charge, LastSteps.discharge, LastSteps.shut]); % Combine all steps
for i = 1:length(StepsOfCycles) % Loop through each step of the cycles
    % Calculate xH2 based on the current state's solution ratio (rs)
    xH2 = rhoGS .* states{i}.rs ./ (rhoOS + rhoGS .* states{i}.rs); 
    plotCellData(model.G, xH2, 'LineStyle', 'none'); 
    time = sum(schedule.step.val(1:i)); 
    title(sprintf('xH2 at Time = %.2f Days', convertTo(time, day)));
    colorbar; 
    pause(0.01);
end
%% Copyright Notice
%
% <html>
% <p><font size="-1">
% Copyright 2009-2025 SINTEF Digital, Mathematics & Cybernetics.
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

