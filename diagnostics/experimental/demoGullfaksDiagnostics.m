clear grdecl
mrstModule add deckformat ad-fi mex mrst-gui coarsegrid diagnostics
prefix = fullfile(ROOTDIR, 'examples', 'data', 'gullfaks', 'E300HIST4');
egrid = readEclipseOutputFileUnFmt([prefix '.EGRID']);
init = readEclipseOutputFileUnFmt([prefix '.INIT']);
cdims = [80 100 52];

grdecl.COORD = egrid.COORD.values;
grdecl.ZCORN = egrid.ZCORN.values;
grdecl.ACTNUM = int32(egrid.ACTNUM.values);
grdecl.cartDims = cdims;
% Create "dumb" grid without splitting disconnected parts
G = processgrid_mex(grdecl);
G = computeGeometry(G);
%% Handle aquifer cell explicitly
active = grdecl.ACTNUM;
active_noaqua = grdecl.ACTNUM;
active_noaqua(280010) = 0;
ind = ismember(find(active), find(active_noaqua));

%% Convert and create rock
grdecl.PERMX = convertFrom(init.PERMX.values(ind), milli*darcy);
grdecl.PERMY = convertFrom(init.PERMY.values(ind), milli*darcy);
grdecl.PERMZ = convertFrom(init.PERMZ.values(ind), milli*darcy);
grdecl.PORO = init.PORO.values(ind);
grdecl.PROPS = [];

grdecl.GRID = grdecl;
rock = initEclipseRock(grdecl);
%% Extract connected subgrid
p = ones(G.cells.num, 1);
p = processPartition(G, p);

connectedCells = find(p == 1);
G = extractSubgrid(G, connectedCells);rock.perm = rock.perm(connectedCells, :);
G.cartDims = cdims;
rock.poro = rock.poro(connectedCells);
G = computeGeometry(G);
% A very small subset is negative because of intersecting pillars. We
% preserve grid ordering and connectivity and simply set it to abs.
G.cells.volumes = abs(G.cells.volumes);
%%
ijk = gridLogicalIndices(G);
W = [];
gc = G.cells.centroids;
[pi,ii] = deal(1);
for i = 1:15:G.cartDims(1)
    for j = 1:20:G.cartDims(2)
        c = ijk{1} == i & ijk{2} == j;
        if any(c)
            c = find(c);
            x = gc(c(1), 1); y = gc(c(1), 2);
            if x > 4.56e5 && x < 4.59e5 && y < 6.79e6
                val = 500*barsa;
                name = ['I' num2str(ii)];
                ii = ii + 1;
            else
                val = 250*barsa;
                name = ['P' num2str(pi)];
                pi = pi + 1;
            end

            W = addWell(W, G, rock, c, 'Type', 'bhp', 'Val', val, 'Name', name);
        end
    end
end
%%
close all
interactiveDiagnostics(G, rock, W);
colormap jet
axis tight off

%%
data.rock = rock;
data.cells = G.cells;

close all
plotToolbar(G, data, 'EdgeColor', 'k', 'EdgeAlpha', .115)
axis tight off
