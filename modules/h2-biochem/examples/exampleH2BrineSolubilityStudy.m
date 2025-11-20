%% MRST Example: Solubility Tables for H2-Brine System with RK-EoS vs SW-EoS vs PR-EoS
%
% This example demonstrates the setup and computation of solubility tables
% for hydrogen (H2) in a brine system using the Redlich-Kwong (RK) Equation
% of State. It also demonstrates the computation of H2 solubility using the
% Peng-Robinson (PR) and the Soreide-Whitson (SW) EoS models with flash calculations,
% with a comparison to tabulated data from ePC-SAFT.
%
%   Ahmed, E., et al. (2024). Phase behavior and black-oil simulations of
%   hydrogen storage in saline aquifers. Advances in Water Resources, 191, 104772.
%   Ahmed, E., et al. (2025). Modeling and simulation of coupled biochemical
%   and two-phase compositional flow in underground hydrogen storage
%
% SEE ALSO:
%   `generateComponentTable`, `generateSolubilityTable
%
%
%{
Copyright 2009-2026 SINTEF Digital, Mathematics & Cybernetics.

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

clear; clc;

mrstModule add biochemistry compositional ad-blackoil ad-core ad-props h2store mrst-gui


% Input Parameters for Temperature, Pressure, and Salinity
min_temp     = 1;                 % [°C]
max_temp     = 99;                % [°C]
min_pressure = 6 * mega * Pascal; % [Pa]
max_pressure = 20 * mega * Pascal;% [Pa]
nbp          = 15;               % Number of pressure points
nbt          = 15;              % Number of temperature points
ms           = 2;                 % Salt molality [mol/kg]
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

%% calculate H2 solubility from PR and SW EoS models (Flash calculations)
%% Define Composition Mixture
% Create a compositional mixture with water and hydrogen components
compFluid = TableCompositionalMixture({'Water', 'Hydrogen'}, {'H2O', 'H2'});
disp(compFluid);

%% Initialize Thermodynamic Model
% Choose an equation of state model for the calculations
eosNamesw = 'sw'; % Soreide-Whitson (SW) model
eosNamepr = 'pr'; % Peng-Robinson (PR) model
eosModelsw = EquationOfStateModel([], compFluid, eosNamesw);
eosModelpr = EquationOfStateModel([], compFluid, eosNamepr);
eosModelsw.msalt=ms;
z0 = [0.8, 0.2]; % Initial composition
pressSWPR=tab_sol.("pressure [Pa]");
tempSWPR=tab_sol.("# temperature [°C]")+ 273.15;

nc = numel(pressSWPR);
namecp = eosModelsw.getComponentNames();
indH2=find(strcmp(namecp,'H2'));
indH2O= find(strcmp(namecp,'H2O'));

[Lsw, xsw, ~] = standaloneFlash(pressSWPR, tempSWPR, z0, eosModelsw);
[Lpr, xpr, ~] = standaloneFlash(pressSWPR, tempSWPR, z0, eosModelpr);
xliqH2sw=xsw(:,indH2);
xliqH2pr=xpr(:,indH2);

%% Also calculate gas phase compositions for SW and PR
% For gas phase H2O content (y_H2O)
% Using the same flash results but extracting vapor phase compositions
% For standaloneFlash, we need to get vapor phase compositions
% Re-run flash to get both phases
[~, xsw_full, ysw_full] = standaloneFlash(pressSWPR, tempSWPR, z0, eosModelsw);
[~, xpr_full, ypr_full] = standaloneFlash(pressSWPR, tempSWPR, z0, eosModelpr);
yliqH2Osw = ysw_full(:, indH2O);
yliqH2Opr = ypr_full(:, indH2O);

%% ========================================================================
% FIGURE 1: H2 Solubility in Brine - All Models Comparison
% ========================================================================
figure('Position', [100, 100, 1200, 800]);
hold on;
plotColors = {'k', 'r', 'b', 'm', 'g'}; % Colors for each pressure level
lineStyles = {'-', '--', '-.', ':', '-'};
kk = 1:int16((n - 1) / 4):n;
legends = cell(20,1); % Initialize legends array
legendIdx = 1;

% Loop over the four selected pressure levels
for j = 1:4
    i = kk(j);
    % Convert the pressure value to MPa and set up legend label
    P_val = unique(pressSWPR(i:n:end));
    P_legend = [num2str(int16(convertTo(P_val, mega * Pascal)))];

    % Plot RK-EoS results with solid lines
    plot(tempSWPR(i:n:n * n), tab_sol.x_H2(i:n:n * n), ...
        'Color', plotColors{j}, 'LineStyle', '-', 'LineWidth', 1.5);
    legends{legendIdx} = ['RK EoS ' P_legend ' MPa'];
    legendIdx = legendIdx + 1;

    % Plot SW-EoS results with dashed lines
    plot(tempSWPR(i:n:n * n), xliqH2sw(i:n:n * n), ...
        'Color', plotColors{j}, 'LineStyle', '--', 'LineWidth', 1.5);
    legends{legendIdx} = ['SW EoS ' P_legend ' MPa'];
    legendIdx = legendIdx + 1;

    % Plot PR-EoS results with dot-dash lines
    plot(tempSWPR(i:n:n * n), xliqH2pr(i:n:n * n), ...
        'Color', plotColors{j}, 'LineStyle', '-.', 'LineWidth', 1.5);
    legends{legendIdx} = ['PR EoS ' P_legend ' MPa'];
    legendIdx = legendIdx + 1;

    % Plot ePC-SAFT results with dotted lines
    plot(state.T(i:n:n * n), state.X_L(i:n:n * n, indexH2), ...
        'Color', plotColors{j}, 'LineStyle', ':', 'LineWidth', 1.5);
    legends{legendIdx} = ['ePC-SAFT ' P_legend ' MPa'];
    legendIdx = legendIdx + 1;
end

xlabel('Temperature (K)', 'FontSize', 14);
ylabel('x_{H2} (mole fraction)', 'FontSize', 14);
legend(legends(1:legendIdx-1), 'FontSize', 9, 'Location', 'best');
title('H2 Solubility in Salt Water: RK-EoS vs. SW-EoS vs. PR-EoS vs. ePC-SAFT', 'FontSize', 12);
grid on;
hold off;

%% ========================================================================
% FIGURE 2: H2O Solubility in Gas Phase - All Models Comparison
% ========================================================================
figure('Position', [100, 100, 1200, 800]);
hold on;
legends = cell(20,1);
legendIdx = 1;

% Loop over the four selected pressure levels
for j = 1:4
    i = kk(j);
    P_val = unique(pressSWPR(i:n:end));
    P_legend = [num2str(int16(convertTo(P_val, mega * Pascal)))];

    % Plot RK-EoS results (from tab_sol)
    plot(tempSWPR(i:n:n * n), tab_sol.y_H2O(i:n:n * n), ...
        'Color', plotColors{j}, 'LineStyle', '-', 'LineWidth', 1.5);
    legends{legendIdx} = ['RK EoS ' P_legend ' MPa'];
    legendIdx = legendIdx + 1;

    % Plot SW-EoS results (from flash)
    plot(tempSWPR(i:n:n * n), yliqH2Osw(i:n:n * n), ...
        'Color', plotColors{j}, 'LineStyle', '--', 'LineWidth', 1.5);
    legends{legendIdx} = ['SW EoS ' P_legend ' MPa'];
    legendIdx = legendIdx + 1;

    % Plot PR-EoS results (from flash)
    plot(tempSWPR(i:n:n * n), yliqH2Opr(i:n:n * n), ...
        'Color', plotColors{j}, 'LineStyle', '-.', 'LineWidth', 1.5);
    legends{legendIdx} = ['PR EoS ' P_legend ' MPa'];
    legendIdx = legendIdx + 1;

    % Plot ePC-SAFT results
    plot(state.T(i:n:n * n), state.X_V(i:n:n * n, indexH2O), ...
        'Color', plotColors{j}, 'LineStyle', ':', 'LineWidth', 1.5);
    legends{legendIdx} = ['ePC-SAFT ' P_legend ' MPa'];
    legendIdx = legendIdx + 1;
end

xlabel('Temperature (K)', 'FontSize', 14);
ylabel('y_{H2O} (mole fraction in gas)', 'FontSize', 14);
legend(legends(1:legendIdx-1), 'FontSize', 9, 'Location', 'best');
title('H2O Vaporization in Hydrogen Gas: RK-EoS vs. SW-EoS vs. PR-EoS vs. ePC-SAFT', 'FontSize', 12);
grid on;
hold off;

%% ========================================================================
% FIGURE 3: Contour Plots - H2 Solubility for All EoS Models
% ========================================================================
% Create a 2x2 subplot for contour comparisons
figure('Position', [100, 100, 1400, 1000]);

% RK-EoS Contour Plot
subplot(2,2,1);
contourf(temperature, convertTo(pressure, mega*Pascal), ...
    reshape(tab_sol.x_H2, [], n), 20, 'LineStyle', 'none');
colorbar;
title('H2 Solubility (RK EoS)', 'FontSize', 12);
xlabel('Temperature (K)', 'FontSize', 10);
ylabel('Pressure (MPa)', 'FontSize', 10);
set(gca, 'FontSize', 10);
grid on;

% SW-EoS Contour Plot
subplot(2,2,2);
contourf(temperature, convertTo(pressure, mega*Pascal), ...
    reshape(xliqH2sw, [], n), 20, 'LineStyle', 'none');
colorbar;
title('H2 Solubility (SW EoS)', 'FontSize', 12);
xlabel('Temperature (K)', 'FontSize', 10);
ylabel('Pressure (MPa)', 'FontSize', 10);
set(gca, 'FontSize', 10);
grid on;

% PR-EoS Contour Plot
subplot(2,2,3);
contourf(temperature, convertTo(pressure, mega*Pascal), ...
    reshape(xliqH2pr, [], n), 20, 'LineStyle', 'none');
colorbar;
title('H2 Solubility (PR EoS)', 'FontSize', 12);
xlabel('Temperature (K)', 'FontSize', 10);
ylabel('Pressure (MPa)', 'FontSize', 10);
set(gca, 'FontSize', 10);
grid on;

% ePC-SAFT Contour Plot - FIXED
subplot(2,2,4);
% Reshape ePC-SAFT H2 solubility data
Z_epcsaft_h2 = reshape(state.X_L(:, indexH2), nbp, nbt);
% Get unique temperature and pressure values
T_unique = state.T(1:n:end);
P_unique = reshape(state.pressure, nbp, nbt);
P_unique_col = P_unique(:, 1);
% Create meshgrid for contour
[T_mesh, P_mesh] = meshgrid(T_unique, P_unique_col);
% Plot contour
contourf(T_mesh, convertTo(P_mesh, mega*Pascal), Z_epcsaft_h2, 20, 'LineStyle', 'none');
colorbar;
title('H2 Solubility (ePC-SAFT)', 'FontSize', 12);
xlabel('Temperature (K)', 'FontSize', 10);
ylabel('Pressure (MPa)', 'FontSize', 10);
set(gca, 'FontSize', 10);
grid on;

sgtitle('H2 Solubility Contour Comparison: All EoS Models', 'FontSize', 14, 'FontWeight', 'bold');

%% ========================================================================
% FIGURE 4: Contour Plots - H2O Vaporization for All EoS Models
% ========================================================================
figure('Position', [100, 100, 1400, 1000]);

% RK-EoS Contour Plot - H2O in Gas
subplot(2,2,1);
contourf(temperature, convertTo(pressure, mega*Pascal), ...
    reshape(tab_sol.y_H2O, [], n), 20, 'LineStyle', 'none');
colorbar;
title('H2O in Gas Phase (RK EoS)', 'FontSize', 12);
xlabel('Temperature (K)', 'FontSize', 10);
ylabel('Pressure (MPa)', 'FontSize', 10);
set(gca, 'FontSize', 10);
grid on;

% SW-EoS Contour Plot - H2O in Gas
subplot(2,2,2);
contourf(temperature, convertTo(pressure, mega*Pascal), ...
    reshape(yliqH2Osw, [], n), 20, 'LineStyle', 'none');
colorbar;
title('H2O in Gas Phase (SW EoS)', 'FontSize', 12);
xlabel('Temperature (K)', 'FontSize', 10);
ylabel('Pressure (MPa)', 'FontSize', 10);
set(gca, 'FontSize', 10);
grid on;

% PR-EoS Contour Plot - H2O in Gas
subplot(2,2,3);
contourf(temperature, convertTo(pressure, mega*Pascal), ...
    reshape(yliqH2Opr, [], n), 20, 'LineStyle', 'none');
colorbar;
title('H2O in Gas Phase (PR EoS)', 'FontSize', 12);
xlabel('Temperature (K)', 'FontSize', 10);
ylabel('Pressure (MPa)', 'FontSize', 10);
set(gca, 'FontSize', 10);
grid on;

% ePC-SAFT Contour Plot - H2O in Gas - FIXED
subplot(2,2,4);
% Reshape ePC-SAFT water vaporization data
Z_epcsaft_h2o = reshape(state.X_V(:, indexH2O), nbp, nbt);
% Use same meshgrid as before
[T_mesh, P_mesh] = meshgrid(T_unique, P_unique_col);
% Plot contour
contourf(T_mesh, convertTo(P_mesh, mega*Pascal), Z_epcsaft_h2o, 20, 'LineStyle', 'none');
colorbar;
title('H2O in Gas Phase (ePC-SAFT)', 'FontSize', 12);
xlabel('Temperature (K)', 'FontSize', 10);
ylabel('Pressure (MPa)', 'FontSize', 10);
set(gca, 'FontSize', 10);
grid on;

sgtitle('H2O Vaporization Contour Comparison: All EoS Models', 'FontSize', 14, 'FontWeight', 'bold');

%% ========================================================================
% FIGURE 5: H2 Density Comparison - All Models
% ========================================================================
% Extract SAFT densities for H2
rhoH2_saft = state.DensH2;
rhoH2_nist = tab_H2.("density [kg/m3]");

% Constants for density calculations
mH2 = 2.016e-3;  % Molar mass of H2 (kg/mol)
Joule = 8.314472;  % Ideal gas constant (J/(mol*K))

% Calculate RK-EoS density from tab_sol
% Using the RK-EoS density from the solubility table
rhoH2_rk = tab_H2.("density [kg/m3]"); % Already computed from NIST for pure H2

% Calculate density using Brill & Beggs correlation
ZH2_corr = calculateBrillBreggsZfactorHydrogen(temperature(:), pressure(:));
rhoH2_corr = pressure(:) .* mH2 ./ (ZH2_corr .* Joule .* temperature(:));
rhoH2_corr = reshape(rhoH2_corr, [], n);

figure('Position', [100, 100, 1200, 800]);
hold on;
legends = cell(16,1);
legendIdx = 1;

for j = 1:4
    i = kk(j);
    P_val = unique(pressSWPR(i:n:end));
    P_legend = [num2str(int16(convertTo(P_val, mega * Pascal)))];

    % RK-EoS (NIST)
    plot(temperature(i:n:n * n), rhoH2_rk(i:n:n * n), ...
        'Color', plotColors{j}, 'LineStyle', '-', 'LineWidth', 1.5);
    legends{legendIdx} = ['RK-EoS/NIST ' P_legend ' MPa'];
    legendIdx = legendIdx + 1;

    % ePC-SAFT
    plot(temperature(i:n:n * n), rhoH2_saft(i:n:n * n, indexH2O), ...
        'Color', plotColors{j}, 'LineStyle', '--', 'LineWidth', 1.5);
    legends{legendIdx} = ['ePC-SAFT ' P_legend ' MPa'];
    legendIdx = legendIdx + 1;

    % Brill & Beggs Correlation
    plot(temperature(i:n:n * n), rhoH2_corr(i:n:n * n), ...
        'Color', plotColors{j}, 'LineStyle', '-.', 'LineWidth', 1.5);
    legends{legendIdx} = ['Brill-Beggs ' P_legend ' MPa'];
    legendIdx = legendIdx + 1;
end

xlabel('Temperature (K)', 'FontSize', 14);
ylabel('H2 Density (kg/m³)', 'FontSize', 14);
legend(legends(1:legendIdx-1), 'FontSize', 9, 'Location', 'best');
title('H2 Density Comparison: RK-EoS/NIST vs. ePC-SAFT vs. Brill-Beggs Correlation', 'FontSize', 12);
grid on;
hold off;

%% ========================================================================
% FIGURE 6: Water Density Comparison - All Models
% ========================================================================
rhoWater_saft = state.DensH2O;
rhoWater_nist = tab_H2O.("density [kg/m3]");
rhoWater_corr = calculateRoweChouWaterDensity(temperature(:), 0);
rhoWater_corr = reshape(rhoWater_corr, [], n);

figure('Position', [100, 100, 1200, 800]);
hold on;
legends = cell(12,1);
legendIdx = 1;

for j = 1:4
    i = kk(j);
    P_val = unique(pressSWPR(i:n:end));
    P_legend = [num2str(int16(convertTo(P_val, mega * Pascal)))];

    % RK-EoS/NIST
    plot(temperature(i:n:n * n), rhoWater_nist(i:n:n * n), ...
        'Color', plotColors{j}, 'LineStyle', '-', 'LineWidth', 1.5);
    legends{legendIdx} = ['RK-EoS/NIST ' P_legend ' MPa'];
    legendIdx = legendIdx + 1;

    % ePC-SAFT
    plot(temperature(i:n:n * n), rhoWater_saft(i:n:n * n), ...
        'Color', plotColors{j}, 'LineStyle', '--', 'LineWidth', 1.5);
    legends{legendIdx} = ['ePC-SAFT ' P_legend ' MPa'];
    legendIdx = legendIdx + 1;
end

% Rowe-Chou Correlation (pressure independent for pure water)
plot(temperature(1:n:n * n), rhoWater_corr(1:n:n * n), ...
    'Color', 'black', 'LineStyle', '-.', 'LineWidth', 2.0);
legends{legendIdx} = 'Rowe-Chou Correlation';
legendIdx = legendIdx + 1;

xlabel('Temperature (K)', 'FontSize', 14);
ylabel('Water Density (kg/m³)', 'FontSize', 14);
legend(legends(1:legendIdx-1), 'FontSize', 9, 'Location', 'best');
title('Water Density Comparison: RK-EoS/NIST vs. ePC-SAFT vs. Rowe-Chou Correlation', 'FontSize', 12);
grid on;
hold off;

%% ========================================================================
% FIGURE 7: Henry's Constant Comparison - All Models (if available)
% ========================================================================
% Calculate Henry's constant from different models
% Henry's constant = fugacity coefficient * pressure / mole fraction

% For RK-EoS, we can compute from tab_sol
H_rk = pressure(:) ./ tab_sol.x_H2(:);
H_rk = reshape(H_rk, [], n);

% For SW-EoS
H_sw = pressure(:) ./ xliqH2sw(:);
H_sw = reshape(H_sw, [], n);

% For PR-EoS
H_pr = pressure(:) ./ xliqH2pr(:);
H_pr = reshape(H_pr, [], n);

% For ePC-SAFT
H_epcsaft = state.pressure(:) ./ state.X_L(:, indexH2);
H_epcsaft = reshape(H_epcsaft, [], n);

figure('Position', [100, 100, 1200, 800]);
hold on;
legends = cell(16,1);
legendIdx = 1;

for j = 1:4
    i = kk(j);
    P_val = unique(pressSWPR(i:n:end));
    P_legend = [num2str(int16(convertTo(P_val, mega * Pascal)))];

    % RK-EoS
    semilogy(temperature(i:n:n * n), H_rk(i:n:n * n) / 1e6, ...
        'Color', plotColors{j}, 'LineStyle', '-', 'LineWidth', 1.5);
    legends{legendIdx} = ['RK EoS ' P_legend ' MPa'];
    legendIdx = legendIdx + 1;

    % SW-EoS
    semilogy(temperature(i:n:n * n), H_sw(i:n:n * n) / 1e6, ...
        'Color', plotColors{j}, 'LineStyle', '--', 'LineWidth', 1.5);
    legends{legendIdx} = ['SW EoS ' P_legend ' MPa'];
    legendIdx = legendIdx + 1;

    % PR-EoS
    semilogy(temperature(i:n:n * n), H_pr(i:n:n * n) / 1e6, ...
        'Color', plotColors{j}, 'LineStyle', '-.', 'LineWidth', 1.5);
    legends{legendIdx} = ['PR EoS ' P_legend ' MPa'];
    legendIdx = legendIdx + 1;

    % ePC-SAFT
    semilogy(state.T(i:n:n * n), H_epcsaft(i:n:n * n) / 1e6, ...
        'Color', plotColors{j}, 'LineStyle', ':', 'LineWidth', 1.5);
    legends{legendIdx} = ['ePC-SAFT ' P_legend ' MPa'];
    legendIdx = legendIdx + 1;
end

xlabel('Temperature (K)', 'FontSize', 14);
ylabel('Henry''s Constant (MPa)', 'FontSize', 14);
legend(legends(1:legendIdx-1), 'FontSize', 9, 'Location', 'best');
title('Henry''s Constant Comparison: All EoS Models', 'FontSize', 12);
grid on;
hold off;

%% ========================================================================
% FIGURE 8: 3D Surface Plot - H2 Solubility Comparison
% ========================================================================
figure('Position', [100, 100, 1400, 1000]);

% Create meshgrid for surface plots - ensure proper dimensions
% temperature is already reshaped as [nbp, nbt] = [15, 15]
% So temperature(1:n:end) gives the first row (all temperatures at min pressure)
T_unique = temperature(1:n:end);  % This gives nbt values (15)
P_unique = convertTo(pressure(1:n:end), mega*Pascal);  % This gives nbp values (15)

% Create proper meshgrid for surfaces
[T_mesh, P_mesh] = meshgrid(T_unique, P_unique);

% RK-EoS Surface
subplot(2,2,1);
% Reshape tab_sol.x_H2 to [nbp, nbt] - it should already be in this shape
Z_rk = reshape(tab_sol.x_H2, nbp, nbt);
surf(T_mesh, P_mesh, Z_rk);
xlabel('Temperature (K)');
ylabel('Pressure (MPa)');
zlabel('x_{H2}');
title('RK EoS');
colormap jet;
colorbar;
view(45, 30);

% SW-EoS Surface
subplot(2,2,2);
Z_sw = reshape(xliqH2sw, nbp, nbt);
surf(T_mesh, P_mesh, Z_sw);
xlabel('Temperature (K)');
ylabel('Pressure (MPa)');
zlabel('x_{H2}');
title('SW EoS');
colormap jet;
colorbar;
view(45, 30);

% PR-EoS Surface
subplot(2,2,3);
Z_pr = reshape(xliqH2pr, nbp, nbt);
surf(T_mesh, P_mesh, Z_pr);
xlabel('Temperature (K)');
ylabel('Pressure (MPa)');
zlabel('x_{H2}');
title('PR EoS');
colormap jet;
colorbar;
view(45, 30);

% ePC-SAFT Surface
subplot(2,2,4);
% Extract ePC-SAFT data and reshape properly
% state.X_L is a vector, need to reshape to [nbp, nbt]
Z_epcsaft = reshape(state.X_L(:, indexH2), nbp, nbt);
% Use the same T_mesh and P_mesh for consistency
surf(T_mesh, P_mesh, Z_epcsaft);
xlabel('Temperature (K)');
ylabel('Pressure (MPa)');
zlabel('x_{H2}');
title('ePC-SAFT');
colormap jet;
colorbar;
view(45, 30);

sgtitle('3D Surface Comparison of H2 Solubility: All EoS Models', 'FontSize', 14, 'FontWeight', 'bold');
%% ========================================================================
% FIGURE 9: Relative Differences Between Models
% ========================================================================
% Calculate relative differences
diff_SW_RK = abs((xliqH2sw - tab_sol.x_H2) ./ tab_sol.x_H2) * 100;
diff_PR_RK = abs((xliqH2pr - tab_sol.x_H2) ./ tab_sol.x_H2) * 100;
diff_SAFT_RK = abs((state.X_L(:, indexH2) - tab_sol.x_H2) ./ tab_sol.x_H2) * 100;

% Reshape for plotting
diff_SW_RK = reshape(diff_SW_RK, [], n);
diff_PR_RK = reshape(diff_PR_RK, [], n);
diff_SAFT_RK = reshape(diff_SAFT_RK, [], n);

figure('Position', [100, 100, 1200, 800]);

% SW vs RK
subplot(1,3,1);
contourf(temperature, convertTo(pressure, mega*Pascal), diff_SW_RK, 20, 'LineStyle', 'none');
colorbar;
title('SW vs RK (\% Difference)', 'FontSize', 12);
xlabel('Temperature (K)', 'FontSize', 10);
ylabel('Pressure (MPa)', 'FontSize', 10);
set(gca, 'FontSize', 10);
grid on;

% PR vs RK
subplot(1,3,2);
contourf(temperature, convertTo(pressure, mega*Pascal), diff_PR_RK, 20, 'LineStyle', 'none');
colorbar;
title('PR vs RK (\% Difference)', 'FontSize', 12);
xlabel('Temperature (K)', 'FontSize', 10);
ylabel('Pressure (MPa)', 'FontSize', 10);
set(gca, 'FontSize', 10);
grid on;

% ePC-SAFT vs RK
subplot(1,3,3);
contourf(temperature, convertTo(pressure, mega*Pascal), diff_SAFT_RK, 20, 'LineStyle', 'none');
colorbar;
title('ePC-SAFT vs RK (\% Difference)', 'FontSize', 12);
xlabel('Temperature (K)', 'FontSize', 10);
ylabel('Pressure (MPa)', 'FontSize', 10);
set(gca, 'FontSize', 10);
grid on;

sgtitle('Relative Differences in H2 Solubility Compared to RK-EoS', 'FontSize', 14, 'FontWeight', 'bold');

%% ========================================================================
% FIGURE 10: Statistical Summary - Box Plot of H2 Solubility
% ========================================================================
figure('Position', [100, 100, 800, 600]);

% Collect all solubility data
all_data = [tab_sol.x_H2(:), xliqH2sw(:), xliqH2pr(:), state.X_L(:, indexH2)];
boxplot(all_data, {'RK-EoS', 'SW-EoS', 'PR-EoS', 'ePC-SAFT'});
ylabel('x_{H2} (mole fraction)', 'FontSize', 12);
title('Statistical Distribution of H2 Solubility Across EoS Models', 'FontSize', 12);
grid on;

%% ========================================================================
% Standard Conditions Data (from original example)
% ========================================================================
% Obtain NIST data for pure components at standard conditions
comp_name = 'H2';
tab_H2_pure = generatePureComponentProperties('min_temp', -3, ...
    'max_temp', 3, 'n_temp', 3, 'comp_name', comp_name, ...
    'outputDisplay', outputDisplay, 'outputPath', outputPath);

comp_name = 'H2O';
tab_H2O_pure = generatePureComponentProperties('min_temp', 15, ...
    'max_temp', 25, 'n_temp', 10, 'comp_name', comp_name, ...
    'outputDisplay', outputDisplay, 'outputPath', outputPath);

% Plot standard conditions data
figure('Position', [100, 100, 1000, 400]);

subplot(1,2,1);
plot(tab_H2O_pure.("# temperature [°C]"), tab_H2O_pure.("density [kg/m3]"), 'b-o', 'LineWidth', 1.5);
xlabel('Temperature [°C]', 'FontSize', 10);
ylabel('Density [kg/m3]', 'FontSize', 10);
title('Water Density at Standard Conditions', 'FontSize', 12);
legend('H2O Density', 'Location', 'Best');
grid on;

subplot(1,2,2);
plot(tab_H2_pure.("# temperature [°C]"), tab_H2_pure.("density [kg/m3]"), 'r-o', 'LineWidth', 1.5);
xlabel('Temperature [°C]', 'FontSize', 10);
ylabel('Density [kg/m3]', 'FontSize', 10);
title('Hydrogen Density at Standard Conditions', 'FontSize', 12);
legend('H2 Density', 'Location', 'Best');
grid on;

sgtitle('Pure Component Properties at Standard Conditions (1 atm)', 'FontSize', 14, 'FontWeight', 'bold');

%% ========================================================================
% Summary Statistics Table
% ========================================================================
fprintf('\n========================================\n');
fprintf('SUMMARY STATISTICS FOR H2 SOLUBILITY\n');
fprintf('========================================\n');
fprintf('Model     | Mean x_H2 | Std Dev   | Min x_H2  | Max x_H2\n');
fprintf('--------------------------------------------------------\n');
fprintf('RK-EoS    | %9.4e | %9.4e | %9.4e | %9.4e\n', ...
    mean(tab_sol.x_H2), std(tab_sol.x_H2), min(tab_sol.x_H2), max(tab_sol.x_H2));
fprintf('SW-EoS    | %9.4e | %9.4e | %9.4e | %9.4e\n', ...
    mean(xliqH2sw), std(xliqH2sw), min(xliqH2sw), max(xliqH2sw));
fprintf('PR-EoS    | %9.4e | %9.4e | %9.4e | %9.4e\n', ...
    mean(xliqH2pr), std(xliqH2pr), min(xliqH2pr), max(xliqH2pr));
fprintf('ePC-SAFT  | %9.4e | %9.4e | %9.4e | %9.4e\n', ...
    mean(state.X_L(:, indexH2)), std(state.X_L(:, indexH2)), ...
    min(state.X_L(:, indexH2)), max(state.X_L(:, indexH2)));
fprintf('========================================\n');

% Force recompute tables
if ~recompute
    disp('You have changed P, T, or ms! Make sure "recompute" is set to true to recalculate tables.');
end