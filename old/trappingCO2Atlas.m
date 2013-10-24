%% Potential CO_2 traps for the CO2 Storage Atlas
% By finding local minima in the surfaces of potential CO2 storage sites,
% we can estimate the locations of structural traps. The grids are created
% from the topographic provided by the Norwegian Petroleum Directorate at
% <http://www.npd.no/en/Publications/Reports/CO2-Storage-Atlas-/>
% The data files are in 
% http://www.npd.no/Global/Engelsk/3-Publications/Reports/CO2-Storage-Atlas/CO2-Atlas-Norwegian-NorthSea-2011.zip
grids = readAtlasGrids();
%% Heightmap of the different potential traps
% We first plot the top surface of all the reservoars we will analyze
% colorized by height of the top surface.
assert(numel(grids)==10)
f = figure( 'Position', [100, 100, 1400, 800]);
for i = 1:numel(grids)
    G = grids{i};
    subplot(2,5,i)
    plotCellData(G, G.cells.z, 'edgec', 'none');
    title(G.name);
    axis tight off
end
%% Find CO2 trapping, and plot for all reservoirs
for i = 1:numel(grids)
    G = grids{i};
    g = G.name;
    % Find trapping
    [z_spill_loc,g_trap,trap_level,z_spill_level,z_spill_loc_level,Gnew]=findTrappingStructure_dfs(G);
    % Estimate trapped volumes
    volumes = estimateStructuralTrapping(G, z_spill_loc);
    %%
    close all;
    f = figure( 'Position', [100, 100, 1400, 800]);
    fprintf('Plotting "%s". %d cells.\n', g, G.cells.num);
    clf;
    colormap jet
    set(f,'name',g,'numbertitle','off')
    plotCellData(G, G.cells.z, 'facea', .3, 'edgea', .05, 'edgec', 'k')
    colorbar;
    plotTrapping(G, z_spill_loc, volumes)
    axis tight
    title(sprintf('Trapping for %s', g));
end
