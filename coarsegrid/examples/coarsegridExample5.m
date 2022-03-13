%% Generate Coarse Grids with Near-Well Refinement
% Pressure gradients and flow rates will typically be much larger near
% wells than inside the reservoir. The accuracy with which we represent the
% flow in and out of wells will to a large extent determine the accuracy of
% an overall simulation and as a result one therefore often desires to have
% higher grid resolution in the near-well zone than inside the reservoir.
%
% In this example, we will use the function 'refineNearWell' to make coarse
% grids with various types of near-well refinement.

mrstModule add coarsegrid incomp

%% Make fine grid and far-field partitioning
% To demonstrate the various coarsening types, we start with a simple
% rectangular domain discretized by a uniform Cartesian grid. As our first
% far-field coarse grid, we use an almost uniform partitioning

G = cartGrid([150 150 6], [100 100 1]*meter);
x = G.nodes.coords(:,1); y = G.nodes.coords(:,2);
G.nodes.coords(:,3) = G.nodes.coords(:,3) - exp( -(x-60).^2./1000 - (y-40).^2/2000); 
clf, plotGrid(G,'EdgeAlpha',.05); view(3);
%%
p0 = partitionUI(G, [5 5 1]);
G = computeGeometry(G);
plotPartition = @(p) plotCellData(G, p, 'EdgeColor','k','EdgeAlpha',.05);

% Make a dummy rock object which is required to set up wells
clear rock
rock.perm = ones(G.cells.num, 1);

% Make two wells: one in the center and one in the SW corner
W = verticalWell([], G, rock, round(G.cartDims(1)/2), ...
                 round(G.cartDims(1)/2), [],'InnerProduct', 'ip_tpf');
W = verticalWell(W, G, rock, 1, 1, [], 'InnerProduct', 'ip_tpf');

% Plot the partitioned grid
clf
plotPartition(p0); outlineCoarseGrid(G, p0)
plotWell(G, W, 'height', 1/2)
axis tight off, view(25, 60), colormap(colorcube(max(p0)))

%% Show two different refinements per well
% The function 'refineNearWell' takes a set of points and partitions these
% according to the distance in the xy-plane from a single well point. Here,
% we will this function to refine the coarse blocks that contains wells.
% For the first well, we will partition the wellblock into five radial and
% six angular sections. The second well lies at the SW corner of the
% wellblock and hence, we only partition it into five radial sections. The
% width of the radial sections is set to decay as log(r).
p = p0;
angSectors = [6 1];
radSectors = [5 5];
for i = 1:numel(W)
    % Find the well blocks
    wc      = W(i).cells(1);
    pt_well = G.cells.centroids(wc, :);

    cells = p == p(wc);
    pts = G.cells.centroids(cells, :);
    out = refineNearWell(pts, pt_well, 'angleBins', angSectors(i), ...
         'radiusBins', radSectors(i), 'logbins', true, 'maxRadius', inf);

    p(cells) = max(p) + out;
end

%%
p1 = partitionCartGrid(G.cartDims,[1 1 3]);
p = p + max(p).*(p1-1);

%%
clf
plotPartition(p); %outlineCoarseGrid(G, p)
CG = generateCoarseGrid(G, p);
plotFaces(CG,1:CG.faces.num,'FaceColor','none','EdgeColor','k','LineWidth',1);
plotWell(G, W, 'height', 1/2)
axis tight off, view(25, 60), colormap(colorcube(max(p)))


%% Increasing number of angular sectors away from the center
% It is also possible to set the number of angular sectors so that it
% increases as we move radially out from the well point
pt_well1 = G.cells.centroids(W(1).cells(1), :);
out = refineNearWell(G.cells.centroids, pt_well1, 'logbins', true, ...
   'angleBins', [2 4 8 16], 'radiusBins', 4, 'maxRadius', inf);

clf
plotPartition(out); outlineCoarseGrid(G, out)
plotWell(G, W(1), 'height', 1/2)
axis tight off, view(25, 60), colormap(colorcube(max(out)))

%% Unstructured coarse grid with radial refinement
% In the last example, we will use METIS to make an unstructured partition
% that adapts to the underlying geology. To this end, we generate a
% lognormal permeability distribution and use the resulting
% transmissibilities as edge-weights in the graph-partitioning algorithm of
% METIS so that it tries to make grid blocks having as homogeneous
% permeability as possible. We then combine this partition with a radial
% refinement. (The result is not guaranteed to give accurate simulations,
% but the example illustrates the flexibility in this type of coarsening.)
rock.perm = logNormLayers(G.cartDims, [300 200]*milli*darcy, 'sz', [51 3 3]);
T = computeTrans(G, rock);
p = partitionMETIS(G, T, 20, 'useLog', true);

p_loc = refineNearWell(G.cells.centroids, pt_well1, 'angleBins', 6, ...
   'radiusBins', 4, 'logbins', true, 'maxRadius', [15 10]*meter);
local = p_loc>0;
p(local) = p_loc(local) + max(p);

clf
plotGrid(G, 'EdgeColor', 'w' , 'EdgeAlpha', .2); outlineCoarseGrid(G, p)
plotWell(G, W(1), 'height', 1/2)
axis tight off, view(25, 60)

axes('Position',[.02 .78 .2 .2]); 
plotCellData(G,log10(rock.perm),'EdgeColor','none'); 
view(25,60), axis tight off; colormap(jet)

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
