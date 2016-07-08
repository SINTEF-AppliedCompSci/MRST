%% Read ECLIPSE output and visualize
% This example reads ECLIPSE unformatted output files from a simulation
% based on the SPE9 benchmark. The MRST grid structure (G) is generated
% using output from *.INIT and *.EGRID. Unconverted grid properties from
% *INIT and *UNRST (restart data) are plotted using plotToolBar. Summary
% data (*UNSMRY) is read and inspected.

mrstModule add mrst-gui

if ~ makeSPE9OutputAvailable,
   error('SPE9Download:Failure', ...
         'Failed to download ECLIPSE output for SPE-9 benchmark case');
end

prefix = fullfile(getDatasetPath('spe9'), 'Simulation-Output', 'SPE9_CP');

%% Read INIT/EGRID-files and construct MRST-grid
init  = readEclipseOutputFileUnFmt([prefix, '.INIT']);
egrid = readEclipseOutputFileUnFmt([prefix, '.EGRID']);
[G, rock] = eclOut2mrst(init, egrid);

%% Plot (unconverted) init data compatible with the grid (G)
figure, plotToolbar(G, init);
axis off vis3d tight, view([1 -1.5 .3]), camproj perspective

%% Read and plot (unconverted) restart data compatible with the grid (G)
% By first reading the restart specification file, we may choose to read
% only portions of the restart data. Here er read all data for for
% time-steps 1 through 20.

spec  = processEclipseRestartSpec(prefix, 'all');
steps = 1:20;
rstrt = readEclipseRestartUnFmt(prefix, spec, steps);

%%
% The data-arrangement in rstrt is not compatible with plotToolbar, hence
% we need to rearrage the data as a struct-array. At the same time we only
% pick data that matches the number of grid cells.
ix = structfun(@(x)numel(x{1}) == G.cells.num, rstrt);
f  = fieldnames(rstrt); f = f(ix);

%%
% Rearrange restart-data to array of structures for plotting
data = cellfun(@(x)rstrt.(x), f, 'UniformOutput', false);
data = cell2struct(vertcat(data{:}), f, 1);

%%
% Plot selected restart data
figure, plotToolbar(G, data);
axis off vis3d tight, view([1, -1.5, 0.3]), camproj perspective

%% Read and inspect (unconverted) summary data
% The summary reader creates a struct smry containg all the data in an
% array accompanied by index-functions for extracting data for a given
% well/property/timestep. The reader requires both the *.UNSMRY and
% *.SMSPEC files.

smry = readEclipseSummaryUnFmt(prefix);

%%
% Display names of all "wells" for which properties are recorded:

smry.WGNAMES

%%
% The field with the expressive name ':+:+:+:+' contains ministep 
% (simulator time-step) times given by the property 'TIME'. 
% Request times of all (':') ministeps: 
time = smry.get(':+:+:+:+', 'TIME', ':');
% Find property names ascossiated with the field 'FIELD'
smry.getKws('FIELD')

%%
% Plot field oil production rate for all time steps
figure, plot(time, smry.get('FIELD', 'FOPR', ':'), 'LineWidth', 2);
xlabel('Time [days]'), ylabel('Field oil production rate [stb/day]');
