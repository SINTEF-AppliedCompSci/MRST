function plotBacterialEffects(scenarios, timeYears, poro0, perm0)
% plotBacterialEffects -- Plot bacterial concentration and clogging effects
%
% SYNOPSIS:
%   plotBacterialEffects(scenarios, timeYears, poro0, perm0)
%
% DESCRIPTION:
%   Generates diagnostic plots for bio-clogging simulations. The function
%   shows the time evolution of:
%       1) Average bacterial concentration
%       2) Porosity reduction (normalized by initial porosity)
%       3) Permeability reduction (normalized by initial permeability)
%
% PARAMETERS:
%   scenarios - Cell array of scenario structs. Each scenario must contain:
%                  * states : cell array of state structs, optionally with
%                             field 'nbact' for bacterial concentration
%                  * model  : MRST model struct (with rock fields poro, perm)
%                  * line   : line style for plotting
%                  * color  : RGB vector for line color
%                  * name   : string for legend
%   timeYears - Vector of times in years for the x-axis.
%   poro0     - Reference porosity values for normalization.
%   perm0     - Reference permeability values for normalization.
%
% RETURNS:
%   none (generates figure).
%
% SEE ALSO:
%   setupBioCloggingModel, plotComponentProfiles

figure('Position', [100, 100, 1000, 800]);

% Identify scenarios with bacterial concentration field
hasBacteria   = cellfun(@(x) isfield(x.states{1}, 'nbact'), scenarios);
bioScenarios  = scenarios(hasBacteria);

if isempty(bioScenarios)
    return;  % Nothing to plot
end

%% 1. Average bacterial concentration
subplot(2, 2, 1);
hold on;
for scenIdx = 1:numel(bioScenarios)
    scen = bioScenarios{scenIdx};
    nSteps = numel(scen.states);
    bactData = zeros(nSteps, 1);

    for step = 1:nSteps
        bactData(step) = mean(scen.states{step}.nbact);
    end

    plot(timeYears, bactData, scen.line, ...
        'Color', scen.color, 'LineWidth', 2, ...
        'DisplayName', scen.name);
end
title('Average Bacterial Concentration');
xlabel('Time (years)');
ylabel('Concentration');
grid on;
legend;

%% 2. Porosity reduction
subplot(2, 2, 2);
hold on;
for scenIdx = 1:numel(bioScenarios)
    scen = bioScenarios{scenIdx};
    nSteps = numel(scen.states);
    poroData = zeros(nSteps, 1);

    for step = 1:nSteps
        currentPoro = scen.model.rock.poro;
        if isa(currentPoro, 'function_handle')
            nbact = scen.states{step}.nbact;
            poroData(step) = mean(currentPoro(1, nbact));
        else
            poroData(step) = mean(currentPoro);
        end
    end

    plot(timeYears, poroData ./ mean(poro0), scen.line, ...
        'Color', scen.color, 'LineWidth', 2, ...
        'DisplayName', scen.name);
end
title('Porosity Reduction (Normalized)');
xlabel('Time (years)');
ylabel('Porosity / Initial Porosity');
grid on;

%% 3. Permeability reduction
subplot(2, 2, 3);
hold on;
for scenIdx = 1:numel(bioScenarios)
    scen = bioScenarios{scenIdx};
    nSteps = numel(scen.states);
    permData = zeros(nSteps, 1);

    for step = 1:nSteps
        currentPerm = scen.model.rock.perm;
        if isa(currentPerm, 'function_handle')
            nbact = scen.states{step}.nbact;
            permData(step) = mean(currentPerm(1, nbact));
        else
            permData(step) = mean(currentPerm(:, 1));
        end
    end

    plot(timeYears, permData ./ mean(perm0), scen.line, ...
        'Color', scen.color, 'LineWidth', 2, ...
        'DisplayName', scen.name);
end
title('Permeability Reduction (Normalized)');
xlabel('Time (years)');
ylabel('Permeability / Initial Permeability');
grid on;

%% Global title
sgtitle('Bio-Clogging Effects');
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