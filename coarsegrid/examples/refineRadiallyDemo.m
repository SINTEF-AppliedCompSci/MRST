%% Refine coarsegrids near wells / specific points in xy-plane
mrstModule add multiscale-devel coarsegrid libgeometry mrst-gui

% Grid
G = cartGrid([150 150 5], [100 100 1]*meter);
p0 = partitionUI(G, [3 3 1]);
G = mcomputeGeometry(G);

% Fake a rock
clear rock
rock.perm = ones(G.cells.num, 1);

% Make two wells
W = [];
W = verticalWell(W, G, rock, round(G.cartDims(1)/2), round(G.cartDims(1)/2), []);
W = verticalWell(W, G, rock, 1, 1, []);

%% Show two different refinements per well
p = p0;
for i = 1:numel(W)
    if i == 1
        sectors = 6;
    else
        sectors = 1;
    end
    wc = W(i).cells(1);

    pt_well = G.cells.centroids(wc, :);

    cells = p == p(wc);
    pts = G.cells.centroids(cells, :);
    out = refineNearWell(pts, pt_well, 'angleBins', sectors, 'radiusBins', 5, 'logbins', true, 'maxRadius', inf);

    p(cells) = max(p) + out;
end
close all
plotToolbar(G, mod(p, 13), 'edgec', 'w' , 'edgea', .2)
outlineCoarseGrid(G, p)
plotWell(G, W, 'height', 1/2)
axis tight off
view(25, 60)
%% Different amount of sectors as we move away from the center
pt_well1 = G.cells.centroids(W(1).cells(1), :);
out = refineNearWell(G.cells.centroids, pt_well1, 'angleBins', [2 4 8 16],...
                                                  'radiusBins', 4,...
                                                  'logbins', true,...
                                                  'maxRadius', inf);

close all
plotToolbar(G, mod(out, 13), 'edgec', 'w' , 'edgea', .2)
outlineCoarseGrid(G, out)
plotWell(G, W(1), 'height', 1/2)
axis tight off
view(25, 60)


%% Use only radius to plot, with a unstructured grid underneath
mrstModule add mrst-experimental
T = computeTrans(G, rock);
p = partitionMETIS(G, T, 10);
%% Refine near the well
p_loc = refineNearWell(G.cells.centroids, pt_well1, 'angleBins', 5, 'radiusBins', 4, 'logbins', true, 'maxRadius', 15*meter);
local = p_loc>0;
p(local) = p_loc(local) + max(p);

close all
plotToolbar(G, mod(p, 13), 'edgec', 'w' , 'edgea', .2)
outlineCoarseGrid(G, p)
plotWell(G, W(1), 'height', 1/2)
axis tight off
view(25, 60)