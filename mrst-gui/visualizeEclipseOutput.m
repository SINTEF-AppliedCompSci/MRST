function [G, data] = visualizeEclipseOutput(prefix)
% Simple, experimental Eclipse output visualization function.
% Intentionally undocumented.

    require deckformat ad-fi
    init = readEclipseOutputFileUnFmt([prefix, '.INIT']);
    grid = readEclipseOutputFileUnFmt([prefix, '.EGRID']);
    [G, rock, N, T] = eclOut2mrst(init, grid);

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
