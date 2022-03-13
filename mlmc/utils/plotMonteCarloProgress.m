function h = plotMonteCarloProgress(out, h, varargin)
%Undocumented Utility Function

%{
Copyright 2009-2022 SINTEF Digital, Mathematics & Cybernetics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}

    if nargin < 2, h = []; end
    opt = struct('color'    , lines(1));
    opt = merge_options(opt, varargin{:});

    if isempty(opt.color)
        opt.color = lines(1);
    end

    if isempty(h)
        df = get(0, 'DefaultFigurePosition');
        h = figure('Position', [df(1:2), 800, 400]);
    else
        set(0, 'CurrentFigure', h); clf(h);
    end

    color = opt.color;
    a = 0.5; pcolor = color.*a + (1-a);
    cost = out.cost.*out.numSamples;
    rmse = out.rmse;

    up = out.estimate + rmse;
    lo = out.estimate - rmse;

    hold on
    ix = 1+isinf(up(1)):numel(cost);
    if isempty(ix), ix = 1; end
    if numel(cost) > 1
        patch([cost(ix); flipud(cost(ix))], [up(ix); flipud(lo(ix))], 1, ...
                               'FaceColor', pcolor, 'EdgeColor', 'none');
        plot(cost(ix), up(ix), 'Color', 'k', 'LineWidth', 1);
        plot(cost(ix), lo(ix), 'Color', 'k', 'LineWidth', 1);
    end

    mrk = {};
    if numel(cost) <= 100, mrk = {'Marker', '.', 'MarkerSize', 25}; end
    plot(cost(ix), out.estimate(ix), 'Color', 'k', 'LineWidth', 2, mrk{:});
    plot([cost(ix(1)), cost(ix(end))], [out.estimate(ix(end)), out.estimate(ix(end))], '--k');
    hold off

    box on; grid on;
    if numel(cost) > 1
        dx = max(cost) - min(cost);
        dy = max(up) - min(lo);
        xlim([min(cost), max(cost)] + [-1,1]*dx*0.01)
        ylim([min(lo)  , max(up)  ] + [-1,1]*dy*0.01)
    end
end
