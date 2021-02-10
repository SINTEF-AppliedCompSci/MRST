% Test marking cut domains from sliceGrid

clear all
close all

N = 5;
G = cartGrid([N, N, N], [20, 20, 20]);
G = computeGeometry(G);
G = addBoundingBoxFields(G);
nCuts = 5;
parentFaces = cell(nCuts, 1);

% Slice randomly
for k = 1:nCuts
    [G, gix] = sliceGrid(G, 10 + 18*rand(1,3)-9, 'normal', 20*rand(1,3)-10,...
        'radius', inf*(1 + 10*rand(1,3)), 'topoSplit', true);
    chk = checkGrid(G);
    assert(chk);

    % Save parent face map
    parentFaces{k} = gix.parent.faces;
end

% Initialize with last face status
faceStatus = gix.new.faces;

% Map face status from first slice to last slice
for m = 2:nCuts
    % Initial faces to map
    fk = find(parentFaces{m-1} == 0);
    for k = m:nCuts
        fk = mapface(parentFaces{k}, fk);
    end
    faceStatus(fk) = 3;
end

% Mark domains
mloop = markCutGrids(G, faceStatus);

% Plot all regions
figure, hold on
N = max(mloop);
color = lines(N);
for k = 1:N
    cells = mloop == k;
    plotCellData(G, mloop, cells, 'facecolor', color(k,:));
end
axis tight, view(3), rotate3d

% Subplot regions
figure
sgtitle 'volume fractions'
V = 20^3;
m = floor(sqrt(N));
n = ceil(N/m);
for k = 1:N
    subplot(m, n, k)
    hold on
    plotGrid(G, 'facecolor', 'none', 'edgealpha', 0.1);
    cells = mloop==k;
    plotGrid(G, cells, 'facecolor', color(k,:));
    title(sprintf('%1.2e', sum(G.cells.volumes(cells))/V));
    axis tight, view(3)
end

function fmap = mapface(parentFaces, f)
% Utility for mapping faces to the parent
    fmap = [];
    for i = 1:numel(f)
        fmap = [fmap; find(parentFaces == f(i))];
    end
end
