%% H2-H2O Vapor-Liquid Equilibrium Calculations
%
% This example calculates the vapor-liquid equilibrium (VLE) for a hydrogen-water
% mixture using two equation of state models: Peng-Robinson (PR) and
% Soreide-Whitson (SW). The solubility of H2 in pure water and NaCl brine
% is computed and compared against experimental data.
%
% The calculations are performed for different pressure, temperature, and
% salinity conditions using standalone flash calculations.
%
% References:
%   Solubility of H2 in water and NaCl brine under subsurface storage conditions
%   https://hal.science/hal-04623907v1, 2023
%
% SEE ALSO:
%   `standaloneFlash`, `EquationOfStateModel`, `TableCompositionalMixture`
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

% Add necessary MRST modules
mrstModule add h2-biochem compositional ad-blackoil ad-core ad-props mrst-gui

%% 1. Define Composition Mixture
% Create a compositional mixture with water and hydrogen components
compFluid = TableCompositionalMixture({'Water', 'Hydrogen'}, {'H2O', 'H2'});
disp('Compositional mixture initialized:');
disp(compFluid);

%% Initialize Thermodynamic Model
% Choose an equation of state model for the calculations
eosModelsw = SoreideWhitsonEos([], compFluid); % Soreide-Whitson (SW) model

eosModelpr = EquationOfStateModel([], compFluid, 'pr'); % Peng Robinson (PR) model

%% 3. Define Test Case Parameters
% Select experimental dataset from literature
% Cases 1-4: H2 solubility in water and NaCl brine at various salinities

z0 = [0.8, 0.2]; % Initial composition [H2O, H2]
caseTest = 1; % Select test case (1, 2, 3, or 4)

switch caseTest
    case 1
        % Pure water system (no salt) at various temperatures
        eosModelsw.msalt = 0;
        eosModelpr.msalt = 0;
        pres = [100.01, 150.01, 200.01, 200.01, 101.11, 101.31, 130.01, ...
            165.01, 199.91, 200.11, 100.01, 100.01, 100.01, 175.11] * barsa;
        Temp = [298.20, 298.05, 298.15, 298.15, 323.55, 323.50, 323.85, 323.55, ...
            323.30, 323.35, 373.85, 373.80, 373.85, 373.65] * Kelvin;
        xliqH2Exp = [0.00135994, 0.00199679, 0.00264397, 0.00263320, 0.00125091, ...
            0.00123925, 0.00159253, 0.00201117, 0.00244186, 0.00245479, ...
            0.00142471, 0.00140332, 0.00140893, 0.00243347];
        caseDescription = 'H2 in pure water at 298-373 K';

    case 2
        % Moderate salinity brine (1 mol/kg) at various temperatures
        eosModelsw.msalt = 1;
        eosModelpr.msalt = 1;
        Temp = [298.20, 298.30, 298.15, 298.30, 323.20, 323.40, 323.35, 323.20, ...
            323.30, 323.30, 323.20, 373.25, 373.40, 373.10, 373.15, 373.45, 373.15, 373.00] * Kelvin;
        pres = [100.71, 150.01, 150.01, 200.01, 100.31, 100.61, 101.01, 150.06, ...
            175.01, 199.91, 200.01, 100.11, 126.01, 150.36, 150.46, 150.71, 175.51, 200.46] * barsa;
        xliqH2Exp = [0.00107012, 0.00159298, 0.00161244, 0.00213575, 0.00102099, ...
            0.00102078, 0.00102426, 0.00151377, 0.00176728, 0.00202590, ...
            0.00204487, 0.00119671, 0.00148492, 0.00175595, 0.00179573, ...
            0.00177350, 0.00204000, 0.00234604];
        caseDescription = 'H2 in brine (1 mol/kg) at 298-373 K';

    case 3
        % High salinity brine (2 mol/kg) at various temperatures
        eosModelsw.msalt = 2;
        eosModelpr.msalt = 2;
        Temp = [298.15, 298.05, 298.05, 323.20, 323.40, 323.35, 323.40, 373.05, 373.20, 373.40] * Kelvin;
        pres = [100.01, 150.01, 200.01, 100.01, 150.01, 150.01, 200.01, 100.01, 150.01, 200.51] * barsa;
        xliqH2Exp = [0.00088640, 0.00132235, 0.00171848, 0.00088260, 0.00131242, ...
            0.00128912, 0.00172402, 0.00099379, 0.00151866, 0.00205031];
        caseDescription = 'H2 in brine (2 mol/kg) at 298-373 K';

    case 4
        % Very high salinity brine (4 mol/kg) at various temperatures
        eosModelsw.msalt = 4;
        eosModelpr.msalt = 4;
        Temp = [298.20, 298.20, 298.15, 323.30, 323.40, 323.40, 373.25, 373.35, 373.15] * Kelvin;
        pres = [100.01, 150.01, 200.01, 100.01, 150.01, 200.01, 100.01, 150.01, 200.01] * barsa;
        xliqH2Exp = [0.00059422, 0.00093595, 0.00121838, 0.00061736, ...
            0.00095752, 0.00129237, 0.00077991, 0.00114509, 0.00157469];
        caseDescription = 'H2 in brine (4 mol/kg) at 298-373 K';
end

fprintf('\n=== Test Case %d: %s ===\n', caseTest, caseDescription);
fprintf('Number of experimental points: %d\n', numel(Temp));
fprintf('Salinity: %.2f mol/kg\n', eosModelsw.msalt);

%% 4. Perform Flash Calculations
% Calculate equilibrium compositions using both EoS models

nc = numel(pres);
namecp = eosModelsw.getComponentNames();
indH2 = find(strcmp(namecp, 'H2'));
indH2O = find(strcmp(namecp, 'H2O'));

fprintf('\nPerforming flash calculations...\n');

% SW-EoS flash calculations
[Lsw, xsw, ysw] = standaloneFlash(pres, Temp, z0, eosModelsw);
xliqH2sw = xsw(:, indH2);

% PR-EoS flash calculations
[Lpr, xpr, ypr] = standaloneFlash(pres, Temp, z0, eosModelpr);
xliqH2pr = xpr(:, indH2);

fprintf('Flash calculations completed.\n');

%% 5. Error Analysis
% Calculate relative errors between model predictions and experimental data

errorSWexp = abs(xliqH2sw - xliqH2Exp') ./ xliqH2Exp';
errorPRexp = abs(xliqH2pr - xliqH2Exp') ./ xliqH2Exp';

fprintf('\n=== Error Analysis ===\n');
fprintf('%-20s %-15s %-15s\n', 'Metric', 'SW-EoS', 'PR-EoS');
fprintf('%-20s %-15s %-15s\n', '------------------', '---------------', '---------------');
fprintf('%-20s %12.8f   %12.8f\n', 'Maximum Error:', max(errorSWexp), max(errorPRexp));
fprintf('%-20s %12.8f   %12.8f\n', 'Minimum Error:', min(errorSWexp), min(errorPRexp));
fprintf('%-20s %12.8f   %12.8f\n', 'Mean Error:', mean(errorSWexp), mean(errorPRexp));
fprintf('%-20s %12.8f   %12.8f\n', 'Median Error:', median(errorSWexp), median(errorPRexp));

%% 6. Group Data by Temperature
% Identify unique temperature groups (round to nearest 0.5 K to handle floating point)
Temp_rounded = round(Temp / 0.5) * 0.5;
uniqueTemps = unique(Temp_rounded);
fprintf('\nFound %d unique temperature groups: ', length(uniqueTemps));
fprintf('%.1f K ', uniqueTemps);
fprintf('\n');

% Create temperature groups
tempGroups = cell(length(uniqueTemps), 1);
presGroups = cell(length(uniqueTemps), 1);
xExpGroups = cell(length(uniqueTemps), 1);
xSwGroups = cell(length(uniqueTemps), 1);
xPrGroups = cell(length(uniqueTemps), 1);

for i = 1:length(uniqueTemps)
    idx = Temp_rounded == uniqueTemps(i);
    tempGroups{i} = Temp(idx);
    presGroups{i} = convertTo(pres(idx), barsa);
    xExpGroups{i} = xliqH2Exp(idx);
    xSwGroups{i} = xliqH2sw(idx);
    xPrGroups{i} = xliqH2pr(idx);

    fprintf('\nTemperature %.1f K: %d experimental points\n', uniqueTemps(i), sum(idx));
end

% Create full pressure array for all points (for main plot)
pres_bar_all = convertTo(pres, barsa);

%% 7. Visualize Results - Multi-Panel Comparison
% Create a 2x2 subplot layout for better understanding
% Define colors for different temperatures
tempColors = lines(length(uniqueTemps));

figure('Position', [100, 100, 1000, 800]);

pres_bar = convertTo(pres, barsa);

% Subplot 1: All data together
subplot(2,2,1);
plot(pres_bar, xliqH2Exp, 'ko', 'MarkerSize', 6, 'DisplayName', 'Exp.');
hold on;
plot(pres_bar, xliqH2sw, 'b-', 'LineWidth', 1.5, 'DisplayName', 'SW');
plot(pres_bar, xliqH2pr, 'r--', 'LineWidth', 1.5, 'DisplayName', 'PR');
xlabel('Pressure (bar)'); ylabel('x_{H2}');
title('All Data');
legend('Location', 'best');
grid on;
hold off;

% Subplot 2: SW-EoS vs Experimental
subplot(2,2,2);
plot(pres_bar, xliqH2Exp, 'ko', 'MarkerSize', 6, 'DisplayName', 'Exp.');
hold on;
plot(pres_bar, xliqH2sw, 'b-', 'LineWidth', 1.5, 'DisplayName', 'SW');
xlabel('Pressure (bar)'); ylabel('x_{H2}');
title(sprintf('SW-EoS Model (Mean Error: %.1f%%)', mean(errorSWexp)*100));
legend('Location', 'best');
grid on;
hold off;

% Subplot 3: PR-EoS vs Experimental
subplot(2,2,3);
plot(pres_bar, xliqH2Exp, 'ko', 'MarkerSize', 6, 'DisplayName', 'Exp.');
hold on;
plot(pres_bar, xliqH2pr, 'r--', 'LineWidth', 1.5, 'DisplayName', 'PR');
xlabel('Pressure (bar)'); ylabel('x_{H2}');
title(sprintf('PR-EoS Model (Mean Error: %.1f%%)', mean(errorPRexp)*100));
legend('Location', 'best');
grid on;
hold off;

% Subplot 4: Error comparison
subplot(2,2,4);
bar([mean(errorSWexp)*100, mean(errorPRexp)*100], 'FaceColor', [0.5 0.5 0.5]);
set(gca, 'XTickLabel', {'SW-EoS', 'PR-EoS'});
ylabel('Mean Error (%)');
title('Model Performance Comparison');
grid on;

sgtitle(sprintf('H_2 Solubility in %s (Salinity: %.1f mol/kg)', caseDescription, eosModelsw.msalt), ...
    'FontSize', 12, 'FontWeight', 'bold');

%% 8. Temperature-Dependent Analysis - Separate Plots for Each Temperature
% Create separate plots for each temperature group

if length(uniqueTemps) > 1
    figure('Position', [100, 100, 1400, 500]);

    for i = 1:length(uniqueTemps)
        subplot(1, length(uniqueTemps), i);

        % Sort by pressure for better visualization
        [pres_sorted, sortIdx] = sort(presGroups{i});
        xExp_sorted = xExpGroups{i}(sortIdx);
        xSw_sorted = xSwGroups{i}(sortIdx);
        xPr_sorted = xPrGroups{i}(sortIdx);

        % Plot experimental data
        plot(pres_sorted, xExp_sorted, 'ko', ...
            'MarkerSize', 8, 'MarkerFaceColor', 'k', 'LineWidth', 1.5, ...
            'DisplayName', 'Experimental');
        hold on;

        % Plot SW-EoS results
        plot(pres_sorted, xSw_sorted, 'b-s', ...
            'LineWidth', 1.5, 'MarkerSize', 6, 'DisplayName', 'SW-EoS');

        % Plot PR-EoS results
        plot(pres_sorted, xPr_sorted, 'r-^', ...
            'LineWidth', 1.5, 'MarkerSize', 6, 'DisplayName', 'PR-EoS');

        xlabel('Pressure (bar)', 'FontSize', 10);
        ylabel('x_{H2}', 'FontSize', 10);
        title(sprintf('T = %.1f K (%d points)', uniqueTemps(i), length(presGroups{i})), ...
            'FontSize', 11, 'FontWeight', 'bold');
        legend('Location', 'best', 'FontSize', 8);
        grid on;
        hold off;
    end

    sgtitle(sprintf('H_2 Solubility by Temperature - %s', caseDescription), ...
        'FontSize', 12, 'FontWeight', 'bold');
end

%% 9. Error Distribution Visualization
% Create histogram of errors

figure('Position', [100, 100, 1000, 400]);

subplot(1,2,1);
histogram(errorSWexp, 10, 'FaceColor', 'b', 'EdgeColor', 'k', 'FaceAlpha', 0.7);
xlabel('Relative Error', 'FontSize', 11);
ylabel('Frequency', 'FontSize', 11);
title(sprintf('SW-EoS Error Distribution\nMean Error: %.2f%%', mean(errorSWexp)*100), 'FontSize', 11);
grid on;

subplot(1,2,2);
histogram(errorPRexp, 10, 'FaceColor', 'r', 'EdgeColor', 'k', 'FaceAlpha', 0.7);
xlabel('Relative Error', 'FontSize', 11);
ylabel('Frequency', 'FontSize', 11);
title(sprintf('PR-EoS Error Distribution\nMean Error: %.2f%%', mean(errorPRexp)*100), 'FontSize', 11);
grid on;

sgtitle(sprintf('Error Distribution for %s', caseDescription), 'FontSize', 12, 'FontWeight', 'bold');

%% 10. Pressure Effect Visualization at Different Temperatures
% Plot pressure effect with different colors for each temperature

if length(uniqueTemps) >= 2
    figure('Position', [100, 100, 1000, 700]);
    hold on;

    for i = 1:length(uniqueTemps)
        % Sort by pressure
        [pres_sorted, sortIdx] = sort(presGroups{i});
        xExp_sorted = xExpGroups{i}(sortIdx);

        % Plot experimental data
        plot(pres_sorted, xExp_sorted, 'o', ...
            'MarkerSize', 6, 'MarkerFaceColor', tempColors(i,:), ...
            'MarkerEdgeColor', 'k', 'LineWidth', 1.5, ...
            'DisplayName', sprintf('Exp. T = %.1f K', uniqueTemps(i)));

        % Plot SW-EoS results with matching color and dashed line
        xSw_sorted = xSwGroups{i}(sortIdx);
        plot(pres_sorted, xSw_sorted, '--', ...
            'Color', tempColors(i,:), 'LineWidth', 1.5, ...
            'DisplayName', sprintf('SW T = %.1f K', uniqueTemps(i)));
    end

    xlabel('Pressure (bar)', 'FontSize', 12);
    ylabel('H_2 Mole Fraction (x_{H2})', 'FontSize', 12);
    title(sprintf('Pressure Effect on H_2 Solubility - %s', caseDescription), ...
        'FontSize', 12, 'FontWeight', 'bold');
    legend('Location', 'best', 'FontSize', 9);
    grid on;
    hold off;
end

%% 11. Summary Statistics Table
fprintf('\n=== Summary Statistics ===\n');
fprintf('%-10s | %-12s | %-10s | %-10s | %-10s\n', 'Model', 'Mean x_H2', 'Std Dev', 'Min x_H2', 'Max x_H2');
fprintf('-----------|--------------|------------|------------|------------\n');
fprintf('%-10s | %12.6f | %10.6f | %10.6f | %10.6f\n', ...
    'SW-EoS', mean(xliqH2sw), std(xliqH2sw), min(xliqH2sw), max(xliqH2sw));
fprintf('%-10s | %12.6f | %10.6f | %10.6f | %10.6f\n', ...
    'PR-EoS', mean(xliqH2pr), std(xliqH2pr), min(xliqH2pr), max(xliqH2pr));
fprintf('%-10s | %12.6f | %10.6f | %10.6f | %10.6f\n', ...
    'Exp.', mean(xliqH2Exp), std(xliqH2Exp), min(xliqH2Exp), max(xliqH2Exp));
fprintf('-----------|--------------|------------|------------|------------\n');

%% 12. Optional: Save Results
% Uncomment to save results to file
% saveResults = true;
% if saveResults
%     outputFile = sprintf('H2_solubility_results_case%d_%s.mat', caseTest, datestr(now, 'yyyymmdd'));
%     save(outputFile, 'pres', 'Temp', 'xliqH2Exp', 'xliqH2sw', 'xliqH2pr', ...
%          'errorSWexp', 'errorPRexp', 'caseDescription');
%     fprintf('\nResults saved to: %s\n', outputFile);
% end

%% End of example
fprintf('\n=== H2-H2O VLE Calculations Completed ===\n');