%% CO2-H2O Vapor-Liquid Equilibrium Calculations
%
% This example calculates the vapor-liquid equilibrium (VLE) for a CO2-water
% mixture using two equation of state models: Peng-Robinson (PR) and
% Soreide-Whitson (SW). The solubility of CO2 in pure water and NaCl brine
% is computed and compared against experimental data.
%
% The calculations are performed for different pressure, temperature, and
% salinity conditions using standalone flash calculations.
%
% References:
%   Chabab, S. (2019). Thermodynamic study of the CO2-H2O system.
%   https://hal.science/hal-02310963v1
%
% SEE ALSO:
%   `standaloneFlash`, `EquationOfStateModel`, `TableCompositionalMixture`
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

clear; clc;

% Add necessary MRST modules
mrstModule add h2-biochem compositional ad-blackoil ad-core ad-props mrst-gui

%% 1. Define Composition Mixture
% Create a compositional mixture with water and carbon dioxide components
compFluid = TableCompositionalMixture({'Water', 'CarbonDioxide'}, {'H2O', 'CO2'});
disp('Compositional mixture initialized:');
disp(compFluid);
%% Initialize Thermodynamic Model
% Choose an equation of state model for the calculations
eosModelsw = SoreideWhitsonEos([], compFluid); % Soreide-Whitson (SW) model

eosModelpr = EquationOfStateModel([], compFluid, 'pr'); % Peng Robinson (PR) model

%% 3. Define Test Case Parameters
% Select experimental dataset from Chabab (2019)
% Cases 1-3: CO2 solubility in NaCl brine at various salinities
% Case 4: CO2 solubility in pure water

z0 = [0.8, 0.2]; % Initial composition [H2O, CO2]
caseTest = 4; % Select test case (1, 2, 3, or 4)

switch caseTest
    case 1
        eosModelsw.msalt=1.13;
        Temp=[323.02, 322.97, 323.03, 323.04, 372.33, 372.31, 372.29, 372.29, 372.25]*Kelvin;
        pres=[53.450, 75.550, 100.350, 145.080, 31.148, 60.500, 108.840, 151.920, 191.980]*barsa;
        xliqCO2Exp=[0.01030, 0.01290, 0.01510, 0.01700,0.00390, 0.00750, 0.01130, 0.01360, 0.01570];


    case 2
        eosModelsw.msalt=1;
        Temp=[373.38, 373.37, 373.41]*Kelvin;
        pres=[16.983, 32.527, 68.182]*barsa;
        xliqCO2Exp=[0.00237, 0.00426, 0.00833];

    case 3
        eosModelsw.msalt=3.01;
        Temp=[342.82, 342.81, 342.82, 372.39, 372.42, 372.41 , 372.43, 372.45, 372.45]*Kelvin;
        pres=[30.391, 72.559, 100.910, 25.556, 71.417, 100.517, 152.433, 199.597, 229.817]*barsa;
        xliqCO2Exp=[0.00441, 0.00880, 0.01057, 0.00292, 0.00707, 0.00878, 0.01141, 0.01258, 0.01337];

    case 4
        eosModelsw.msalt=0;
        Temp=[323.2, 323.2, 323.2, 323.2, 323.2, 323.2, 333.2, 333.2, 333.2,...
            333.2, 333.2, 333.2, 353.1, 353.1, 353.1, 353.1, 353.1, 353.1]*Kelvin;
        pres=[40.5, 60.6, 80.8, 100.9, 121, 141.1, 40.5, 60.6, 80.8, 100.9,...
            121, 141.1, 40.5, 60.6, 80.8, 100.9, 121, 131]*barsa;
        xliqCO2Exp=[0.0109, 0.0161, 0.019, 0.0205, 0.0214, 0.0217, 0.0096,...
            0.0138, 0.0166, 0.0186, 0.0201, 0.0208, 0.008, 0.0114, 0.014, ...
            0.016, 0.0176, 0.0184];

end

fprintf('\n=== Test Case %d: %s ===\n', caseTest, caseDescription);
fprintf('Number of experimental points: %d\n', numel(Temp));
fprintf('Salinity: %.2f mol/kg\n', eosModelsw.msalt);

%% 4. Perform Flash Calculations
% Calculate equilibrium compositions using both EoS models

nc = numel(pres);
namecp = eosModelsw.getComponentNames();
indCO2 = find(strcmp(namecp, 'CO2'));
indH2O = find(strcmp(namecp, 'H2O'));

fprintf('\nPerforming flash calculations...\n');

% SW-EoS flash calculations
[Lsw, xsw, ysw] = standaloneFlash(pres, Temp, z0, eosModelsw);
xliqCO2sw = xsw(:, indCO2);

% PR-EoS flash calculations
[Lpr, xpr, ypr] = standaloneFlash(pres, Temp, z0, eosModelpr);
xliqCO2pr = xpr(:, indCO2);

fprintf('Flash calculations completed.\n');

%% 5. Error Analysis
% Calculate relative errors between model predictions and experimental data

errorSWexp = abs(xliqCO2sw - xliqCO2Exp') ./ xliqCO2Exp';
errorPRexp = abs(xliqCO2pr - xliqCO2Exp') ./ xliqCO2Exp';

fprintf('\n=== Error Analysis ===\n');
fprintf('%-20s %-15s %-15s\n', 'Metric', 'SW-EoS', 'PR-EoS');
fprintf('%-20s %-15s %-15s\n', '------------------', '---------------', '---------------');
fprintf('%-20s %12.8f   %12.8f\n', 'Maximum Error:', max(errorSWexp), max(errorPRexp));
fprintf('%-20s %12.8f   %12.8f\n', 'Minimum Error:', min(errorSWexp), min(errorPRexp));
fprintf('%-20s %12.8f   %12.8f\n', 'Mean Error:', mean(errorSWexp), mean(errorPRexp));
fprintf('%-20s %12.8f   %12.8f\n', 'Median Error:', median(errorSWexp), median(errorPRexp));

%% 6. Visualize Results - Main Comparison Plot
% Create comparison plot for CO2 solubility

figure('Position', [100, 100, 900, 700]);
hold on;

% Convert pressure to bar for plotting
pres_bar = convertTo(pres, barsa);

% Plot experimental data with color coding by temperature
for i = 1:length(pres)
    plot(pres_bar(i), xliqCO2Exp(i), 'o', ...
        'MarkerSize', 8, 'MarkerFaceColor', markerColors(i,:), ...
        'MarkerEdgeColor', 'k', 'LineWidth', 1);
end

% Plot SW-EoS results
plot(pres_bar, xliqCO2sw, 'b-s', 'LineWidth', 1.5, 'MarkerSize', 6, ...
    'DisplayName', 'SW-EoS Model');

% Plot PR-EoS results
plot(pres_bar, xliqCO2pr, 'r-^', 'LineWidth', 1.5, 'MarkerSize', 6, ...
    'DisplayName', 'PR-EoS Model');

xlabel('Pressure (bar)', 'FontSize', 12);
ylabel('CO_2 Mole Fraction in Liquid Phase (x_{CO2})', 'FontSize', 12);
title(sprintf('CO_2 Solubility in %s', caseDescription), 'FontSize', 12, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 10);
grid on;

% Add salinity information as text annotation
annotation('textbox', [0.15, 0.85, 0.3, 0.05], ...
    'String', sprintf('Salinity: %.2f mol/kg', eosModelsw.msalt), ...
    'FontSize', 10, 'BackgroundColor', 'white', 'EdgeColor', 'none');

% Add temperature legend if multiple temperatures
uniqueTemps = unique(Temp);
if length(uniqueTemps) > 1
    % Create custom legend for temperatures
    tempColors = lines(length(uniqueTemps));
    for i = 1:length(uniqueTemps)
        plot(nan, nan, 'o', 'MarkerSize', 8, 'MarkerFaceColor', tempColors(i,:), ...
            'MarkerEdgeColor', 'k', 'DisplayName', sprintf('Exp. T = %.1f K', uniqueTemps(i)));
    end
    legend('Location', 'best', 'FontSize', 9);
end

hold off;

%% 7. Error Distribution Visualization
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

%% 8. Temperature-Dependent Analysis (if multiple temperatures)
% Create separate plots for each temperature

if length(uniqueTemps) > 1
    figure('Position', [100, 100, 1200, 400]);

    for i = 1:min(3, length(uniqueTemps))
        subplot(1, min(3, length(uniqueTemps)), i);
        tempIdx = Temp == uniqueTemps(i);

        % Convert pressure to bar for this subset
        pres_subset_bar = convertTo(pres(tempIdx), barsa);

        % Plot experimental data for this temperature
        plot(pres_subset_bar, xliqCO2Exp(tempIdx), 'ko', ...
            'MarkerSize', 8, 'MarkerFaceColor', 'k', 'LineWidth', 1.5, ...
            'DisplayName', 'Experimental');
        hold on;

        % Plot SW-EoS results
        plot(pres_subset_bar, xliqCO2sw(tempIdx), 'b-s', ...
            'LineWidth', 1.5, 'MarkerSize', 6, 'DisplayName', 'SW-EoS');

        % Plot PR-EoS results
        plot(pres_subset_bar, xliqCO2pr(tempIdx), 'r-^', ...
            'LineWidth', 1.5, 'MarkerSize', 6, 'DisplayName', 'PR-EoS');

        xlabel('Pressure (bar)', 'FontSize', 10);
        ylabel('x_{CO2}', 'FontSize', 10);
        title(sprintf('T = %.1f K', uniqueTemps(i)), 'FontSize', 11, 'FontWeight', 'bold');
        legend('Location', 'best', 'FontSize', 8);
        grid on;
        hold off;
    end

    sgtitle(sprintf('CO_2 Solubility by Temperature - %s', caseDescription), ...
        'FontSize', 12, 'FontWeight', 'bold');
end

%% 9. Summary Statistics Table
fprintf('\n=== Summary Statistics ===\n');
fprintf('%-10s | %-12s | %-10s | %-10s | %-10s\n', 'Model', 'Mean x_CO2', 'Std Dev', 'Min x_CO2', 'Max x_CO2');
fprintf('-----------|--------------|------------|------------|------------\n');
fprintf('%-10s | %12.6f | %10.6f | %10.6f | %10.6f\n', ...
    'SW-EoS', mean(xliqCO2sw), std(xliqCO2sw), min(xliqCO2sw), max(xliqCO2sw));
fprintf('%-10s | %12.6f | %10.6f | %10.6f | %10.6f\n', ...
    'PR-EoS', mean(xliqCO2pr), std(xliqCO2pr), min(xliqCO2pr), max(xliqCO2pr));
fprintf('%-10s | %12.6f | %10.6f | %10.6f | %10.6f\n', ...
    'Exp.', mean(xliqCO2Exp), std(xliqCO2Exp), min(xliqCO2Exp), max(xliqCO2Exp));
fprintf('-----------|--------------|------------|------------|------------\n');

%% 10. Optional: Save Results
% Uncomment to save results to file
% saveResults = true;
% if saveResults
%     outputFile = sprintf('CO2_solubility_results_case%d_%s.mat', caseTest, datestr(now, 'yyyymmdd'));
%     save(outputFile, 'pres', 'Temp', 'xliqCO2Exp', 'xliqCO2sw', 'xliqCO2pr', ...
%          'errorSWexp', 'errorPRexp', 'caseDescription');
%     fprintf('\nResults saved to: %s\n', outputFile);
% end

%% End of example
fprintf('\n=== CO2-H2O VLE Calculations Completed ===\n');