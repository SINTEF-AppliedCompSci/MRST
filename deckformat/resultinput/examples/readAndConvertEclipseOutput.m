% interactiveDiagnosticsFromEclipseOutput
mrstModule add mrst-gui ad-core
prefix =  '/Users/steink/data/opm-data/spe9/eclipse-simulation/SPE9_CP';

% Read INIT/EGRID-files and construct MRST-grid
init = readEclipseOutputFileUnFmt([prefix, '.INIT']);
egrid = readEclipseOutputFileUnFmt([prefix, '.EGRID']);
[G, rock] = eclOut2mrst(init, egrid);
G = computeGeometry(G);

% Convert restart s
states = convertRestartToStates(prefix, G);


% Fetch wellsols and restart-step times
wellSols = cellfun(@(x)x.wellSol, states, 'UniformOutput', false);
time     = cellfun(@(x)x.time, states);

figure, plotToolbar(G, states);
plotWellData(G, wellSols{1}, 'labelFontSize', 10, 'color', [.3 .3 .3]);
axis off vis3d tight, view([1 -.5 1]), camproj perspective

plotWellSols(wellSols, time);

