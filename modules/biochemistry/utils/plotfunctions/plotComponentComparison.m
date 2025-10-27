function plotComponentComparison(scenarios, timeYears, idxH2, idxCO2, idxCH4)
% plotComponentComparison -- Compare H2 loss, CO2 consumption, and CH4 production
%
% SYNOPSIS:
%   plotComponentComparison(scenarios, timeYears, idxH2, idxCO2, idxCH4)
%
% DESCRIPTION:
%   Generates a three-panel figure comparing component mass changes relative
%   to the abiotic scenario. The plots show:
%       1) H2 loss (%)
%       2) CO2 consumption (%)
%       3) CH4 production (%)
%
% PARAMETERS:
%   scenarios - Cell array of scenario structs. Each scenario must include:
%                  * states : cell array of state structs with field
%                             FlowProps.ComponentTotalMass
%                  * line   : line style for plotting
%                  * color  : RGB vector for line color
%                  * name   : string for legend
%   timeYears - Vector of times in years for the x-axis.
%   idxH2     - Index of hydrogen component in ComponentTotalMass.
%   idxCO2    - Index of carbon dioxide component in ComponentTotalMass.
%   idxCH4    - Index of methane component in ComponentTotalMass.
%
% RETURNS:
%   none (generates figure).
%
% SEE ALSO:
%   plotComponentProfiles, plotBacterialEffects

% Find abiotic scenario for reference
scenarioNames = cellfun(@(x) x.name, scenarios, 'UniformOutput', false);
abioticIdx = find(strcmp(scenarioNames, 'Abiotic'), 1);
if isempty(abioticIdx)
    warning('No abiotic scenario found; cannot compute comparisons.');
    return;
end

abioticStates = scenarios{abioticIdx}.states;
nSteps = numel(abioticStates);

% Compute total component masses for abiotic scenario using cellfun
totalH2_abiotic  = cellfun(@(s) sum(s.FlowProps.ComponentTotalMass{idxH2}), abioticStates);
totalCO2_abiotic = cellfun(@(s) sum(s.FlowProps.ComponentTotalMass{idxCO2}), abioticStates);
totalCH4_abiotic = cellfun(@(s) sum(s.FlowProps.ComponentTotalMass{idxCH4}), abioticStates);

%% Create figure
figure('Position', [100, 100, 1200, 400]);

%% 1. H2 Loss
subplot(1, 3, 1); hold on;
for scenIdx = setdiff(1:numel(scenarios), abioticIdx)
    scen = scenarios{scenIdx};
    totalH2 = cellfun(@(s) sum(s.FlowProps.ComponentTotalMass{idxH2}), scen.states);
    H2_loss = ((totalH2_abiotic - totalH2) ./ totalH2_abiotic) * 100;
    plot(timeYears, H2_loss, scen.line, ...
        'Color', scen.color, 'LineWidth', 2, 'DisplayName', scen.name);
end
title('H2 Loss (%)'); xlabel('Time (years)'); ylabel('H2 Loss (%)'); grid on; legend;

%% 2. CO2 Consumption
subplot(1, 3, 2); hold on;
for scenIdx = setdiff(1:numel(scenarios), abioticIdx)
    scen = scenarios{scenIdx};
    totalCO2 = cellfun(@(s) sum(s.FlowProps.ComponentTotalMass{idxCO2}), scen.states);
    CO2_consumption = ((totalCO2_abiotic - totalCO2) ./ totalCO2_abiotic) * 100;
    plot(timeYears, CO2_consumption, scen.line, ...
        'Color', scen.color, 'LineWidth', 2, 'DisplayName', scen.name);
end
title('CO2 Consumption (%)'); xlabel('Time (years)'); ylabel('CO2 Consumption (%)'); grid on;

%% 3. CH4 Production
subplot(1, 3, 3); hold on;
for scenIdx = setdiff(1:numel(scenarios), abioticIdx)
    scen = scenarios{scenIdx};
    totalCH4 = cellfun(@(s) sum(s.FlowProps.ComponentTotalMass{idxCH4}), scen.states);
    CH4_production = ((totalCH4 - totalCH4_abiotic) ./ totalCH4_abiotic) * 100;
    plot(timeYears, CH4_production, scen.line, ...
        'Color', scen.color, 'LineWidth', 2, 'DisplayName', scen.name);
end
title('CH4 Production (%)'); xlabel('Time (years)'); ylabel('CH4 Production (%)'); grid on;

%% Global title
sgtitle('Component Mass Comparison Relative to Abiotic Scenario');
end

%{
Copyright 2009-2025 SINTEF Digital, Mathematics & Cybernetics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify it under the terms of the GNU 
General Public License as published by the Free Software Foundation, either version 3 of 
the License, or (at your option) any later version.

MRST is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without 
even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the 
GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with MRST. If not, 
see <http://www.gnu.org/licenses/>.
%}