function fn = getPlotAfterStepCO2(state0, model, az, el)
% Get a function that allows for dynamic plotting in `simulateScheduleAD`
% for the assessment of CO2 in the ad-micp module. The parameters az and
% el are the azimuth and elevation angles for the current axes.
%
% This function is modified from a file in the MATLAB Reservoir Simulation
% Toolbox (MRST), see
%   mrst/modules/solvent/utils/getPlotAfterStepSolvent.m
%
% We refer to that function for a complete commented version of the file.
% In this file, we comment on some of the lines.

%{
Partial copyright 2009-2021, SINTEF Digital, Mathematics & Cybernetics.
Partial copyright 2021-2026, NORCE Research AS, Computational Geosciences and
Modeling.

This file is part of the ad-micp module.

ad-micp is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

ad-micp is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this file. If not, see <http://www.gnu.org/licenses/>.
%}
    grid = model.G;
    defaultFigurePosition = get(0, 'defaultFigurePosition');

    resultFigure = figure( ...
        'Position', [defaultFigurePosition(1 : 2), 800, 1200]);

    colorMap = flipud(jet);
    numberOfColors = size(colorMap, 1);

    co2Axes = axes('Parent', resultFigure);
    co2Plot = plotCellData( ...
        grid, state0.s(:, 2), 'EdgeColor', 'none');

    title(co2Axes, 'CO$_2$ [-]', 'Interpreter', 'latex');
    colorbar(co2Axes);
    view(co2Axes, az, el);

    firstColor = round(70 * numberOfColors / 256);
    lastColor = round(100 * numberOfColors / 256);

    colormap(co2Axes, colorMap(firstColor : lastColor, :));
    caxis(co2Axes, [0 1]);

    fn = @(model, states, reports, solver, schedule, simtime) ...
        afterStepFunction( ...
        model, states, reports, solver, schedule, simtime, ...
        resultFigure, co2Plot);
end

function [model, states, reports, solver, ok] = afterStepFunction( ...
        model, states, reports, solver, schedule, simtime, ...
        resultFigure, co2Plot)

    computedStates = cellfun( ...
        @(state) ~isempty(state), states);

    computedReports = cellfun( ...
        @(report) ~isempty(report), reports);

    currentStateIndex = find( ...
        computedStates, 1, 'last');

    if isempty(currentStateIndex)
        ok = true;
        return
    end

    currentState = states{currentStateIndex};
    completedSimulationTime = simtime(computedReports);

    [~, boundaryCells] = boundaryFaces(model.G);

    if ishandle(resultFigure) && ishandle(co2Plot)
        set(0, 'CurrentFigure', resultFigure);

        co2Plot.CData = ...
            currentState.s(boundaryCells, 2);

        drawnow;
    end

    ok = simulationRuntimePanel( ...
        model, states, reports, solver, schedule, ...
        completedSimulationTime);
end
