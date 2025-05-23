%% MRST Example: Solubility Tables for H2-Brine System
%
% This example demonstrates the setup and computation of solubility tables
% for hydrogen (H2) in a brine system using the Redlich-Kwong (RK) Equation
% of State, with a comparison to tabulated data from ePC-SAFT, and
% Henry-Setschnow correlation as presented in:
%
%   Ahmed, E., et al. (2024). Phase behavior and black-oil simulations of
%   hydrogen storage in saline aquifers. Advances in Water Resources, 191, 104772.
%
% SEE ALSO:
%   `generateComponentTable`, `generateSolubilityTable
%
%
%{
Copyright 2009-2025 SINTEF Digital, Mathematics & Cybernetics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software Foundation,
either version 3 of the License, or (at your option) any later version.

MRST is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST. If not, see <http://www.gnu.org/licenses/>.
%}

% Input Parameters for Temperature, Pressure, and Salinity
min_temp     = 1;                 % [°C]
max_temp     = 99;                % [°C]
min_pressure = 6 * mega * Pascal; % [Pa]
max_pressure = 20 * mega * Pascal;% [Pa]
nbp          = 15;                % Number of pressure points
nbt          = 15;                % Number of temperature points
ms           = 0;                 % Salt molality [mol/kg]
outputDisplay = true;            % Set to true to display generated tables
recompute = true;                 % recompuete solubility  tables

%% Notice on Computational Cost
warning('ComputationalCost:Medium', ...
    ['Please be advised that for large nbp and nbt this example often takes a long time ', ...
    'to run: this script will extract data from https://webbook.nist.gov']);

% Define target output directory and create if it doesn't exist
outputPath = fullfile(mrstOutputDirectory(), 'UHS_PVT', 'H2SolubilityTable');

% This configuration prepares solubility data for the H2-brine system
% under RK-EOS in the context of saline aquifer storage.

% Generate H2O pure Component Table from NIST
comp_name = 'H2O';
disp(['Generating component table for: ', comp_name]);
tab_H2O = generateComponentProperties('min_temp',min_temp, 'max_temp',max_temp, ...
    'n_temp', nbt, 'min_press',min_pressure, 'max_press', max_pressure, ...
    'n_press',nbp, 'comp_name', comp_name,'outputDisplay', outputDisplay, ...
    'outputPath',outputPath);

% Generate H2 pure Component Table from NIST
pause(0.1);  % Ensure smooth execution between commands
comp_name = 'H2';
disp(['Generating component table for: ', comp_name]);
tab_H2 = generateComponentProperties('min_temp',min_temp, 'max_temp',max_temp, ...
    'n_temp', nbt, 'min_press',min_pressure, 'max_press', max_pressure, ...
    'n_press',nbp, 'comp_name', comp_name,'outputDisplay', outputDisplay, ...
    'outputPath',outputPath);

% We use the Redlich Kwong Eos to obtain the solubility
disp('Generating solubility table for H2-brine mixture...');

%[tab_sol, status_sol, file_path_sol] = generateSolubilityTable(min_temp, ... 
% max_temp, min_pressure, max_pressure, nbp, nbt, ms, outputPath);

tab_sol= generateH2WaterSolubilityTable('min_temp',min_temp, 'max_temp',max_temp, ...
    'n_temp', nbt,'min_press',min_pressure, 'max_press',max_pressure, ...
    'n_press', nbp, 'ms', ms,'outputDisplay', outputDisplay,'outputPath',outputPath, ...
    'reCompute', recompute);

% Load ePC-SAFT Data
epcsaft = load(fullfile(ROOTDIR,'..','modules','h2store','examples','data',...
    'PcSaftSolubilityTable','ePcSaftH2BrineData.mat'));
state = epcsaft.state;

% Define constants
n = nbt;                         % Number of temperature points
indexH2 = 2;                      % Index of H2 component in liquid phase
indexH2O = 1;                   % Index of water in vapor phase

pressure = reshape(tab_sol.("pressure [Pa]"), [], n);
temperature = reshape(tab_sol.("# temperature [°C]"), [], n) + 273.15;  % Convert to Kelvin

%% Solubility of H2 in Brine as a Function of Pressure and Temperature

figure;
hold on;
plotColors = {'b', 'r', 'g', 'm'}; % Colors for each pressure level
kk = 1:int16((n - 1) / 4):n;
legends = cell(8,1); % Initialize legends array
% Loop over the four selected pressure levels
for j = 1:4
    i = kk(j);
    % Convert the pressure value to MPa and set up legend label
    P_val = unique(state.pressure(i:n:end));
    P_legend = [num2str(int16(convertTo(P_val, mega * Pascal)))];
    % Plot ePC-SAFT results with solid lines
    plot(state.T(i:n:n * n), state.X_L(i:n:n * n, indexH2), ...
        'Color', plotColors{j}, 'LineStyle', '-', 'LineWidth', 1.5);
    legends{2*j-1} = ['ePC-SAFT ' P_legend ' MPa']; % Add legend entry
    % Plot RK-EoS results with dashed lines for the same pressures
    plot(state.T(i:n:n * n), tab_sol.x_H2(i:n:n * n), ...
        'Color', plotColors{j}, 'LineStyle', '--', 'LineWidth', 1.5);
    legends{2*j} = ['RK ' P_legend ' MPa']; % Add legend entry
end

xlabel('Temperature (K)', 'FontSize', 14);
ylabel('x_{H2}', 'FontSize', 14);
legend(legends, 'FontSize', 10, 'Location', 'best');
title('H2 Solubility in Brine: ePC-SAFT vs. RK-EoS', 'FontSize', 12);
hold off;

%% Solubility of H2O in Hydrogen as a Function of Pressure and Temperature

figure;
hold on;
plotColors = {'b', 'r', 'g', 'm'}; % Colors for each pressure level
kk = 1:int16((n - 1) / 4):n;
legends = cell(8,1); % Initialize legends array
% Loop over the four selected pressure levels
for j = 1:4
    i = kk(j);
    % Convert the pressure value to MPa and set up legend label
    P_val = unique(state.pressure(i:n:end));
    P_legend = [num2str(int16(convertTo(P_val, mega * Pascal)))];
    % Plot ePC-SAFT results with solid lines
    plot(state.T(i:n:n * n), state.X_V(i:n:n * n, indexH2O), ...
        'Color', plotColors{j}, 'LineStyle', '-', 'LineWidth', 1.5);
    legends{2*j-1} = ['ePC-SAFT ' P_legend ' MPa']; % Add legend entry
    % Plot RK-EoS results with dashed lines for the same pressures
    plot(state.T(i:n:n * n), tab_sol.y_H2O(i:n:n * n), ...
        'Color', plotColors{j}, 'LineStyle', '--', 'LineWidth', 1.5);
    legends{2*j} = ['RK ' P_legend ' MPa']; % Add legend entry
end

xlabel('Temperature (K)', 'FontSize', 14);
ylabel('y_{H2O}', 'FontSize', 14);
legend(legends, 'FontSize', 10, 'Location', 'best');
title('H2O Solubility in Gas: ePC-SAFT vs. RK-EoS', 'FontSize', 12);
hold off;

%% Plot the Amount of H2 in Water
figure;
contourf(temperature, convertTo(pressure, mega*Pascal), ...
    reshape(tab_sol.x_H2, [], n), 'LineStyle', 'none');
colorbar; % Add a color bar for reference
title('H2 Solubility in Water (RK eos)', 'FontSize', 14);
xlabel('Temperature (K)', 'FontSize', 14);
ylabel('Pressure (MPa)', 'FontSize', 14);
set(gca, 'FontSize', 12); % Set font size for axes
grid on; % Add grid for better readability

%% Plot the Amount of H2O in H2
figure;
contourf(temperature, convertTo(pressure, mega*Pascal), ...
    reshape(tab_sol.y_H2O, [], n), 'LineStyle', 'none');
colorbar; % Add a color bar for reference
title('H2O Vaporization in H2 (RK eos)', 'FontSize', 14);
xlabel('Temperature (K)', 'FontSize', 14);
ylabel('Pressure (MPa)', 'FontSize', 14);
set(gca, 'FontSize', 12); % Set font size for axes
grid on; % Add grid for better readability

%% We can also use Henry-Setschnow correlation
sol_H2 = HenrySetschenowH2BrineEos(temperature, 0, pressure);
figure;
contourf(temperature, convertTo(pressure, mega*Pascal), ...
    reshape(sol_H2.x_H2, [], n), 'LineStyle', 'none');
colorbar; % Add a color bar for reference
title('H2 Solubility in Water (HS eos)', 'FontSize', 14);
xlabel('Temperature (K)', 'FontSize', 14);
ylabel('Pressure (MPa)', 'FontSize', 14);
set(gca, 'FontSize', 12); % Set font size for axes
grid on; % Add grid for better readability

%% Density Comparison for Pure Components: H2 and H2O
% This section compares the densities of pure hydrogen (H2) and water (H2O)
% obtained from various sources: ePC-SAFT calculations, NIST data, and the
% Rowe-Chou correlation.

% Extract SAFT densities for H2 and H2O from state variables
rhoH2_saft = state.DensH2;
rhoWater_saft = state.DensH2O;
rhoWater_nist = tab_H2O.("density [kg/m3]");
rhoH2_nist = tab_H2.("density [kg/m3]");

% Constants
mH2 = 2.016e-3;  % Molar mass of H2 (kg/mol)
mWater = 0.0180153;  % Molar mass of H2O (kg/mol)
Joule = 8.314472;  % Ideal gas constant (J/(mol*K))

% Extract RK-EoS table columns and convert units for comparison
Temperature_K_ = tab_sol.("# temperature [°C]") + 273.15;  % Convert to Kelvin
Pressure_Pa_ = tab_sol.("pressure [Pa]");

% Calculate density of hydrogen using the Brill and Beggs correlation
ZH2_corr = calculateBrillBreggsZfactorHydrogen(Temperature_K_, Pressure_Pa_);
rhoH2_corr = Pressure_Pa_ .* mH2 ./ (ZH2_corr .* Joule .* Temperature_K_);

% Calculate density of pure water using the Rowe-Chou correlation
rhoWater_corr = calculateRoweChouWaterDensity(Temperature_K_, 0);

%% Plotting
figure;
hold on;
% Loop over the four selected pressure levels
legends = cell(12,1);
for j = 1:4
    i = kk(j);
    % Convert the pressure value to MPa and set up legend label
    P_val = unique(state.pressure(i:n:end));
    P_legend = [num2str(int16(convertTo(P_val, mega * Pascal)))];
    
    % Plot ePC-SAFT results for H2 with solid lines
    plot(Temperature_K_(i:n:n * n), rhoH2_saft(i:n:n * n, indexH2O), ...
        'Color', plotColors{j}, 'LineStyle', '-', 'LineWidth', 1.5);
    legends{3*j-2} = ['ePC-SAFT ' P_legend ' MPa']; % Add legend entry    
    
    % Plot NIST results for H2 with dashed lines for the same pressures
    plot(Temperature_K_(i:n:n * n), rhoH2_nist(i:n:n * n), ...
        'Color', plotColors{j}, 'LineStyle', '--', 'LineWidth', 1.5);
    legends{3*j-1} = ['NIST ' P_legend ' MPa']; % Add legend entry        
    
    % Plot Brill & Beggs correlation results for H2 with dash-dot lines
    plot(Temperature_K_(i:n:n * n), rhoH2_corr(i:n:n * n), ...
        'Color', plotColors{j}, 'LineStyle', '-.', 'LineWidth', 2.0);
    legends{3*j} = ['Corr. ' P_legend ' MPa']; % Add legend entry       

end

xlabel('Temperature (K)', 'FontSize', 14);
ylabel('Density (kg/m³)', 'FontSize', 14);
legend(legends, 'FontSize', 10, 'Location', 'best');
title('H2 Density', 'FontSize', 12);
hold off;

%% Plotting Water Density
figure;
hold on;
legends = cell(9,1);
% Loop over the four selected pressure levels
for j = 1:4
    i = kk(j);
    % Convert the pressure value to MPa and set up legend label
    P_val = unique(state.pressure(i:n:end));
    P_legend = [num2str(int16(convertTo(P_val, mega * Pascal)))];   
    
    % Plot ePC-SAFT results for water density with solid lines
    plot(Temperature_K_(i:n:n * n), rhoWater_saft(i:n:n * n), ...
        'Color', plotColors{j}, 'LineStyle', '-', 'LineWidth', 1.5);
    legends{2*j-1} = ['ePC-SAFT ' P_legend ' MPa']; % Add legend entry      
    
    
    % Plot NIST results for water density with dashed lines for the same pressures
    plot(Temperature_K_(i:n:n * n), rhoWater_nist(i:n:n * n), ...
        'Color', plotColors{j}, 'LineStyle', '--', 'LineWidth', 1.5);
    legends{2*j} = ['NIST ' P_legend ' MPa']; % Add legend entry         

end

% Add plot of the Rowe-Chou correlation which depend only on (T,salt)
plot(Temperature_K_(1:n:n * n), rhoWater_corr(1:n:n * n), ...
    'Color', 'black', 'LineStyle', '-.', 'LineWidth', 2.0);
RoweCorr = 'Corr. (Rowe-Chou)';
legends{9} = RoweCorr;
xlabel('Temperature (K)', 'FontSize', 14);
ylabel('Density (kg/m³)', 'FontSize', 14);
legend(legends, 'FontSize', 10, 'Location', 'best');
title('Water Density', 'FontSize', 12);
hold off;

%% Obtain NIST data for pure components at std conditions
% We follow ECLIPSE for the std conditions
% Retrieve NIST data for pure hydrogen (H2)
% In Eclipse std conditions for gases are 1 atm and 0 C
comp_name = 'H2';
tab_H2_pure =  generatePureComponentProperties('min_temp',-3, ...
    'max_temp',3, 'n_temp', 3, 'comp_name', comp_name, ...
    'outputDisplay', outputDisplay,'outputPath',outputPath);

% Retrieve NIST data for pure water (H2O)
comp_name = 'H2O';
% In Eclipse std conditions for liquids are 1 atm and 15 C to 25 C
tab_H2O_pure =  generatePureComponentProperties('min_temp',15, ...
    'max_temp',25, 'n_temp', 10, 'comp_name', comp_name, ...
    'outputDisplay', outputDisplay,'outputPath',outputPath);

% Note: The data is calculated at standard atmospheric pressure (1 atm)
figure;
title('H2O density at std conditions', 'FontSize', 12);
% Plotting H2O density
plot(tab_H2O_pure.("# temperature [°C]"), tab_H2O_pure.("density [kg/m3]"));
xlabel('Temperature [°C]', 'FontSize', 10);
ylabel('Density [kg/m3]', 'FontSize', 10);
legend('H2O Density', 'Location', 'Best');
grid on;
figure;

% Plotting H2 density
title('H2 density at std conditions', 'FontSize', 12);
plot(tab_H2_pure.("# temperature [°C]"), tab_H2_pure.("density [kg/m3]"));
xlabel('Temperature [°C]', 'FontSize', 10);
ylabel('Density [kg/m3]', 'FontSize', 10);
legend('H2 Density', 'Location', 'Best');
grid on;

% Force recompute tables
if ~recompute
    disp('You have changed P, T, or ms! Make sure "recompute" is set to true to recalculate tables.');
end