%% Faulted 2.5D reservoir grid
% In this example, we use the upr module to construct a 2.5 D faulted
% reservoir grid.
mrstModule add upr

%% Plotting functionality
fig2D = @() figure('Position', [0,0,800,500]);
fig3D = @() figure('Position', [0,0,1300, 650]);
alpha = 0.6;
cmap  = jet*alpha + (1-alpha);
setAxProps = @(ax) set(ax, 'View'              , [-10, 45]         , ...
                           'PlotBoxAspectRatio', [4.40, 1.86, 1.00], ...
                           'Projection'        , 'Perspective'     , ...
                           'Box'               , 'on'              , ...
                           'ColorMap'          , cmap              );

%% Load constraining curves
% We start by loading a data structure that contains points describing the
% outline, faults, and well positions
pth = fullfile(mrstPath('upr'), 'datasets');
rp  = load(fullfile(pth, 'reservoirPoints.mat'));
rp  = rp.reservoirPoints;

fig2D(), hold on
plot(rp.outline(:,1), rp.outline(:,2), 'k')
plotLinePath(rp.faultLines,'b');
plotLinePath(rp.wellLines,'.r', 'markerSize', 20);
box on, axis equal tight

%% Generate 2D PEBI grid
% We construct a 2D PEBI grid from the points using pebiGrid, with
% refinement around the wells
rng(2019)
n   = 25; % Approximate number of cells in x-direction
L   = max(rp.outline);
G2D = pebiGrid2D(max(L)/n, L          , ...
    'polybdr'        , rp.outline   , ... % Outline
    'faceConstraints', rp.faultLines, ... % Fault lines
    'FCFactor'       , 0.8          , ... % Relative size of fault cells
    'cellConstraints', rp.wellLines , ... % Well coordinates
    'CCRefinement'   , true         , ... % Refine
    'CCFactor'       , 0.1          , ... % Relative size of well cells
    'interpolateFC'  , true         , ... % Interpolate along fault lines
    'CCEps'          , 0.08*max(L)  );    % Refinement transition
G2D = removeShortEdges(G2D, 1); % The grid may contain very short edges.
                                % We remove these
                                
%% Plot 2D grid
% Well cells are identified with G.cells.tag = true. We will need the well
% number later, so we find and store these.
wellNo2D                = nan(G2D.cells.num,1);
wellNo2D(G2D.cells.tag) = 1:nnz(G2D.cells.tag);
fig2D(), plotGrid(G2D); axis equal tight, box on                 % Grid
plotFaces(G2D, G2D.faces.tag, 'edgeColor', 'b', 'lineWidth', 2); % Faults
plotGrid(G2D, G2D.cells.tag , 'faceColor', 'r');                 % Wells

%% Find regions
% The faults divide the reservoir into six distinct regions. We identify
% these using functionality from the coarsegrid module
mrstModule add coarsegrid
p = ones(G2D.cells.num,1);
p = processPartition(G2D, p, find(G2D.faces.tag));

%% Make 2.5D reservoir model
% We construct a volumetric reservoir model by extruding the 2D grid using
% makeLayeredGrid
nLayers = 23; % Number of layers
layerThickness = ones(nLayers,1)*5 + (rand(nLayers,1)*2 - 1)*2;
G0           = makeLayeredGrid(G2D, layerThickness);
G0           = computeGeometry(G0);
G0.cells.tag = repmat(G2D.cells.tag, nLayers, 1);
G0.faces.tag = false(G0.faces.num,1);
G0.faces.tag(abs(G0.faces.normals(:,3))<.01) = repmat(G2D.faces.tag, nLayers, 1);
wellNo       = repmat(wellNo2D, nLayers, 1);
layerID      = reshape(repmat(1:nLayers, G2D.cells.num, 1), [], 1);
compartID    = repmat(p,nLayers,1);

%% Plot the resulting layered grid
prm = randperm(nLayers)';
fig3D(); plotCellData(G0, prm(layerID), 'edgealpha', 0.2)
outlineCoarseGrid(G0, compartID,'EdgeColor','w','LineWidth',2);
setAxProps(gca), camlight();
colormap(jet)
axis off

%% Remove cells
% To add more realism to the reservoir, we mimic erosion and inaccessible
% parts of the formation by removing cells
bnds = [0 8; 1 7; 1 7; 5 3; 2 6; 5 3];
flag = false(G0.cells.num,1);
for i=1:6
    flag = flag | ((compartID==i) & ...
        (layerID<=bnds(i,1) | layerID>=(nLayers-bnds(i,2))));
end
[G, cellMap] = removeCells(G0, flag);
G      = computeGeometry(G);
layer  = layerID(cellMap);
parts  = compartID(cellMap);
wellNo = wellNo(cellMap);

%% Populate with petrophysical properties
% We populate the grid with petrophysical properties drawn from a layered,
% lognormal, isotropic permeability field. To fake faults with
% displacement, the permeability inside each compartment is sampled from
% the same cube. This produces a plausible layering structure, but does not
% preserve the areal correlation within a single geological layer on
% opposite sides of a fault. 
permMean = [10, 912, 790, 90, 10];   % Mean permeability in each layer
N        = [90, 30, nLayers];        % Num of points in each axial direcion
ind      = [1,5,13,15,20,nLayers+1]; % Layer indices
K        = reshape(logNormLayers(N, permMean, 'indices', ind), N);

perm  = nan(G.cells.num,1);
for i = 1:6
    idx = parts==i;
    perm(idx) = sampleFromBox(G, K, find(idx))*milli*darcy;
end
lperm = log10(perm);
poro = (lperm - min(lperm))./(max(lperm) - min(lperm)).*0.7 + 0.1;
rock = makeRock(G, perm, poro);

%% Plot permeability
fig3D(); plotCellData(G, log10(rock.perm), 'edgealpha', 0.2)
setAxProps(gca), camlight();

%% Shift vertical coordinates
% Finally, we shift the vertical coordinates of the grid to mimc geological
% activity.
x    = G.nodes.coords;
xmax = max(x);
xmin = min(x);
xr   = (x - mean(x))./((xmax - xmin)/2)*2; % Transform to reference domain
z    = peaks(xr(:,1), xr(:,2));            % Vertical shift
x(:,3) = x(:,3) - mean(x(:,3));
% We shift the vertical coordinates by using the built-in MATLAB function
% peaks, and also make it thicker along the y axis
x(:,3) = (x(:,3) + z*7).*(x(:,2)./xmax(2)+1).^2*0.25;
x(:,3) = x(:,3) - min(x(:,3)) + 1000*meter; % Normal reservoir depth
G.nodes.coords = x;
G = computeGeometry(G);

%% Make wells
W = [];
for wNo = 1:max(wellNo)
    cells = find(wellNo == wNo);
    W     = addWell(W, G, rock, cells, 'name', '');
end
pw = @() plotWell(G, W, 'color', 'k', 'height', 25, 'LineWidth', 3);

%% Plot the final result
% View 1: with wells
fig3D(), plotCellData(G, log10(rock.perm), 'edgealpha', 0.2); pw();
setAxProps(gca), camlight, view([-10,45]);

% View 2: exploded view of compartments
fig3D()
xmid = mean(G.cells.centroids);
for i = 1:6
    g = G;
    xlmid = mean(G.cells.centroids(parts==i,:));
    v = xlmid - xmid;
    g.nodes.coords = g.nodes.coords + 0.3*v;
    plotCellData(g, log10(rock.perm), parts==i, 'edgealpha', 0.2);
end
setAxProps(gca); camlight
axis tight
