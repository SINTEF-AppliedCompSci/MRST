function plotPressureAndSaturation(scenarios, timeYears)
% plotPressureAndSaturation -- Plot average pressure and gas saturation
%
% SYNOPSIS:
%   plotPressureAndSaturation(scenarios, timeYears)
%
% DESCRIPTION:
%   Generates a two-panel figure showing:
%       1) Average reservoir pressure (in bar)
%       2) Average gas saturation
%   for all provided scenarios over time. Each scenario is overlaid for comparison.
%
% PARAMETERS:
%   scenarios - Cell array of scenario structs. Each scenario must include:
%                  * states : cell array of state structs with fields
%                             'pressure' and 's' (saturation)
%                  * line   : line style for plotting
%                  * color  : RGB vector for line color
%                  * name   : string for legend
%   timeYears - Vector of times in years for the x-axis.
%
% RETURNS:
%   none (generates figure).
%
% SEE ALSO:
%   plotComponentProfiles, plotBacterialEffects

figure('Position', [100, 100, 1200, 500]);

%% 1. Average Pressure
subplot(1, 2, 1);
hold on;
for scenIdx = 1:numel(scenarios)
    scen = scenarios{scenIdx};
    nSteps = numel(scen.states);
    pressureData = zeros(nSteps, 1);

    for step = 1:nSteps
        pressureData(step) = mean(scen.states{step}.pressure) / barsa;
    end

    plot(timeYears, pressureData, scen.line, ...
        'Color', scen.color, 'LineWidth', 2, ...
        'DisplayName', scen.name);
end
title('Average Pressure');
xlabel('Time (years)');
ylabel('Pressure (bar)');
grid on;
legend;

%% 2. Average Gas Saturation
subplot(1, 2, 2);
hold on;
for scenIdx = 1:numel(scenarios)
    scen = scenarios{scenIdx};
    nSteps = numel(scen.states);
    satData = zeros(nSteps, 1);

    for step = 1:nSteps
        satData(step) = mean(scen.states{step}.s(:, 2)); % Gas saturation
    end

    plot(timeYears, satData, scen.line, ...
        'Color', scen.color, 'LineWidth', 2, ...
        'DisplayName', scen.name);
end
title('Average Gas Saturation');
xlabel('Time (years)');
ylabel('Gas Saturation');
grid on;

%% Global title
sgtitle('Pressure and Saturation Profiles');
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