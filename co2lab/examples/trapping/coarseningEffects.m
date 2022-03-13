%% Johansen formation: the effect of coarsening on trapping
% When working with subsea reservoirs coarsening will always be a factor:
%
% * Simulations on very large grids become prohibitively expensive in terms
% of computing power very fast.
% * All grids are essentially coarse approximations of highly complex
% geometry. Different grids have a different scale for each cell depending
% on what kind of data was used to produce the grid in the first place and
% the intended usage of the datasets.
%
% This example demonstrates the effects of geometry coarsening on a model
% from the CO2 Storage Atlas grids, with a special focus on structural
% trapping. To this end, we generate six realizations of the Johansen
% formation: The first is the full dataset, the second coarsened by a
% factor 2 in both i and j directions, the third a factor 3, and so on.
% This gives a set of grids where the finest has approximately 80,000
% cells, while the most coarse has ~2,000 cells. In the finest grid, each
% cell has a resolution of 500x500 m^2 per cell, while the coarsest has
% 3000x3000 m^2. The resolution in both grids is fairly large in compared
% with typical simulation grids in petroleum application.

mrstModule add co2lab coarsegrid deckformat

N = 6;
[Grids, res] = deal(cell(numel(N),1));
for i = 1:N
    fprintf(1,'\nLoading Johansen formation (coarsening factor: %d)...\n', i);
    gr = getAtlasGrid('Johansenfm', 'coarsening', i);
    G = processGRDECL(gr{1});
    % Create top surface grid
    Gt = topSurfaceGrid(computeGeometry(G(1)));
    Grids{i} = Gt;
    % Create trapping and store volumes of each trap
    res{i} = trapAnalysis(Gt, true);
    res{i}.volumes = volumesOfTraps(Gt, res{i}, []);
end

%% Plot a subset of the formation with different degree of coarsening
% We first define a subdomain consisting of a minimum and maximum value for
% both x and y coordinates which is then plotted on the fine grid.
G = Grids{1};

subdomain = [5.25e5, 6.70e6; 5.30e5, 6.75e6];

x = G.cells.centroids(:,1);
y = G.cells.centroids(:,2);

subset = x > subdomain(1,1) & x < subdomain(2,1) &...
         y > subdomain(1,2) & y < subdomain(2,2);
clf;
plotCellData(G, G.cells.z,'EdgeColor','none')
plotGrid(G, subset, 'EdgeColor', 'None', 'FaceColor', 'black', 'FaceAlpha', .5)

axis tight off

%%
% We can then loop over the grids while plotting the parts of the grid
% within the bounding box with a small height offset to visualize the fine
% scale details that are lost when coarsening the fine scale data. Note
% that there are several features in the topmost finest grid that have
% disappeared when coarsening: In the coarsest plot at the bottom almost
% all details are lost except for the rightmost fault.
clf
colors = {'b', 'r', 'g', 'c', 'm', 'y'};
colorize = @(i) colors{mod(i, length(colors)) + 1};
for i = 1:N

   G = Grids{i};
   G.nodes.z = G.nodes.z + 1000*i;
    
   x = G.cells.centroids(:,1);
   y = G.cells.centroids(:,2);
   subset = x > subdomain(1,1) & x < subdomain(2,1) &...
            y > subdomain(1,2) & y < subdomain(2,2);
         
   plotGrid(G, subset, 'facec', colorize(i), 'edgec', 'k', 'edgea', .3)
   fprintf('Grid with coarsening %d has %d fine cells and z-standard deviation %2.4f\n', ...
      i, G.cells.num, std(G.cells.z))
end
title('Coarsening of subset')
legend(regexp(num2str((1:N)*500), '  ', 'split'), 'Location', 'EastOutside')
view(86, 12)
axis tight off

%% The effect of coarsening on trapping analysis
% It is obvious that fine structural details are lost when coarsening the
% grids. The coarsening operation acts as a smoother on the grid, removing
% ridges, folds and oscillations that are present on a shorter wavelength
% than the coarse cells. Unfortunately, these small oscillations are
% especially interesting for CO2 migration studies: Small local height
% maxima can divert small "rivers" of CO2 and act as structural traps.
%
% To demonstrate that these traps are lost when coarsening, we plot the
% structural traps estimated by trapAnalysis for the different grids. Note
% that several smaller traps are removed as the coarsening increases, which
% can be shown statistically by noting that the mean of the trap volume
% quickly increases as the smaller traps are smoothed away. 
%
% The total trapping volume also changes as the coarsening is increased: In
% the beginning the volume increases as the largest traps become slightly
% larger due to the lower resolution. Later on, the total volume shrinks as
% smaller traps are removed entirely.

% Plot the traps
clf
defaultpos = get(0, 'DefaultFigurePosition');
set(gcf, 'Position', [defaultpos(1:2) - [0 100], 1000 500]);
subplot('position', [.025 .025 .65 .925]);
for i = 1:N
    G = Grids{i};
    G.nodes.z = G.nodes.z + 1000*i;
    
    x = G.cells.centroids(:,1);
    y = G.cells.centroids(:,2);
    subset = x > subdomain(1,1) & x < subdomain(2,1) &...
             y > subdomain(1,2) & y < subdomain(2,2);
         
    plotGrid(G, subset, 'facec', colorize(i), 'facea', .2, ...
       'edgec', colorize(i), 'edgea', .9)

    tr = res{i};
    G_flat = flattenTraps(G, tr);
    G_flat.nodes.z = G_flat.nodes.z + 1000*i;
    plotGrid(G_flat, subset & tr.traps ~= 0, ...
       'FaceColor', colorize(i), 'EdgeColor', 'none')
end
% view(-30, 60)
view(86, 12)
axis tight off
title('Traps for successively coarsed grids')

% Plot the total volume as subplot
axes('position', [.75 .55 .225 .375]);
hold on
vol = cellfun(@(x) sum(x.volumes), res);
for i = 1:N
    bar(i, vol(i), colorize(i))
end
title('   Total trap volume')
set(gca, 'Color', get(gcf, 'Color'), ...
   'XTickLabel', regexp(num2str((1:N)*500), '  ', 'split'))
axis tight

% Plot the average volume as subplot
axes('Position', [0.75 .05 .225 .375])
hold on
mvol = cellfun(@(x) mean(x.volumes), res);
for i = 1:N
    bar(i, mvol(i), colorize(i))
end
title('     Average trap volume')
set(gca, 'Color', get(gcf, 'Color'), ...
   'XTickLabel', regexp(num2str((1:N)*500), '  ', 'split'))
axis tight

%%
% <html>
% <p><font size="-1">
% Copyright 2009-2022 SINTEF Digital, Mathematics & Cybernetics.
% </font></p>
% <p><font size="-1">
% This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).
% </font></p>
% <p><font size="-1">
% MRST is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% </font></p>
% <p><font size="-1">
% MRST is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% </font></p>
% <p><font size="-1">
% You should have received a copy of the GNU General Public License
% along with MRST.  If not, see
% <a href="http://www.gnu.org/licenses/">http://www.gnu.org/licenses</a>.
% </font></p>
% </html>
