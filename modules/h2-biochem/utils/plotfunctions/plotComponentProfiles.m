function plotComponentProfiles(scenarios, componentNames, timeYears)
% plotComponentProfiles -- Plot mole fraction evolution for all scenarios
%
% SYNOPSIS:
%   plotComponentProfiles(scenarios, componentNames, timeYears)
%
% DESCRIPTION:
%   Produces a multi-panel plot showing the time evolution of the average
%   mole fraction for each component across all provided scenarios. Each
%   subplot corresponds to one component, and different scenarios are
%   overlaid for comparison.
%
% PARAMETERS:
%   scenarios      - Cell array of scenario structs. Each scenario must
%                    include fields:
%                       * states : cell array of state structs with field
%                                  `components`
%                       * line   : line style for plotting
%                       * color  : RGB vector for line color
%                       * name   : string for legend
%   componentNames - Cell array of component names (strings).
%   timeYears      - Vector of times in years for the x-axis.
%
% RETURNS:
%   none (generates figure).
%
% SEE ALSO:
%   plotWellSols, plotToolbar

% Create figure
figure('Position', [100, 100, 1200, 800]);
nComp = numel(componentNames);

% Loop over all components
for compIdx = 1:nComp
    subplot(2, 2, compIdx);
    hold on;

    % Loop over all scenarios
    for scenIdx = 1:numel(scenarios)
        scen   = scenarios{scenIdx};
        nSteps = numel(timeYears);
        compData = zeros(nSteps, 1);

        % Extract average mole fraction for this component
        for step = 1:nSteps
            comps = scen.states{step}.components;
            if iscell(comps)
                % Components stored as cell array
                compData(step) = mean([comps{:, compIdx}]);
            else
                % Components stored as matrix
                compData(step) = mean(comps(:, compIdx));
            end
        end

        % Plot scenario data
        plot(timeYears, compData, scen.line, ...
            'Color', scen.color, ...
            'LineWidth', 2, ...
            'DisplayName', scen.name);
    end

    % Format subplot
    title(componentNames{compIdx});
    xlabel('Time (years)');
    ylabel('Average Mole Fraction');
    grid on;
    legend('Location', 'best');
end

% Add global title
sgtitle('Component Mole Fraction Evolution');
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