function [G, data] = visualizeEclipseOutput(prefix)
% Simple, experimental Eclipse output visualization function.
% Intentionally undocumented.

%{
Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.

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

    require deckformat

    init = readEclipseOutputFileUnFmt([prefix, '.INIT']);
    grid = readEclipseOutputFileUnFmt([prefix, '.EGRID']);
    [G, rock, N, T] = initGridFromEclipseOutput(init, grid);

    try
        [rstrt, rsspec] = readRestartLocal(prefix);
    catch e
        try
            [rstrt, rsspec] = readEclipseRestartFmt(prefix);
        catch e
            e.message
            return
        end
    end

    G = computeGeometry(G);
    data = restartToStructArray(rstrt, rock);
    figure; plotToolbar(G, data)
    axis tight off
end

function data = restartToStructArray(rstrt, rock)
    data = [];
    for i = 1:min(structfun(@numel, rstrt))
        d = structfun(@(x) x{i}, rstrt, 'UniformOutput', false);
        d.rock = rock;
        if isfield(rstrt, 'SWAT') && isfield(rstrt, 'SGAS')
            d.saturations = [1 - d.SWAT - d.SGAS, d.SGAS, d.SWAT];
        end
        data = [data; d];
    end
end
