% readAndVisualizeEclipseOutput

mrstModule add mrst-gui

prefix =  '/Users/steink/data/opm-data/spe9/eclipse-simulation/SPE9_CP';

% Read INIT/EGRID-files and construct MRST-grid
init = readEclipseOutputFileUnFmt([prefix, '.INIT']);
egrid = readEclipseOutputFileUnFmt([prefix, '.EGRID']);
[G, rock] = eclOut2mrst(init, egrid);

% Plot init-data corresponing to G
figure, plotToolbar(G, init);
axis off vis3d tight, view([1 -1.5 .3]), camproj perspective

% Read
spec = processEclipseRestartSpec(prefix, 'all');
steps = 1:20;
rstrt = readEclipseRestartUnFmt(prefix, spec, steps);

% pick only fields matching to G.cells.num
ix = structfun(@(x)numel(x{1})==G.cells.num, rstrt);
f  = fieldnames(rstrt); f = f(ix);
% rearrange restart-data to array of structs for plotting
data = cellfun(@(x)rstrt.(x), f, 'UniformOutput', false);
data = cell2struct(vertcat(data{:}), f, 1);

% Plot restart (unconverted) data
figure, plotToolbar(G, data);
axis off vis3d tight, view([1 -1.5 .3]), camproj perspective

% Read summary
smry = readEclipseSummaryUnFmt(prefix);
% display "wells"
smry.WGNAMES
% Field ':+:+:+:+' contains ministep time-data. Request times of all
% ministeps
time = smry.get(':+:+:+:+', 'TIME', ':');
% Find data recorded for 'FIELD'
smry.getKws('FIELD')
% Plot field oil production rate for all time steps
figure, plot(time, smry.get('FIELD', 'FOPR', ':'), 'LineWidth', 2);
xlabel('Time [days]'), ylabel('Field oil production rate [stb/day]');


 