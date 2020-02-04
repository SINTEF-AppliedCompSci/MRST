%% Gridding a subset of SPE10 model 2 using the UPR module
% In this example, we look at how the upr module can be used to grid a
% subset of one layer of SPE10 model 2.

%% Add modules
mrstModule add upr spe10 mrst-gui;

%% Plotting helper
spe10Fig = @(name) figure('Position', [0, 0, 600, 600], 'name', name);

%% Get reference model
% We extract the southern part of layer 49 from SPE10 model 2, which
% comprises a complex pattern of fluvial sandstone channels on a background
% of low-permable mudrock.
[~, modelRef, ~] = setupSPE10_AD('layers', 49, 'J', 1:110);
G     = modelRef.G;
rock  = modelRef.rock;
fluid = modelRef.fluid;

%% Get permeability
% We will use the permeability to construct a grid that honors the
% high permeability contrasts
lperm = log10(modelRef.rock.perm(:,1));
lperm = lperm - min(lperm);
% Plot the base 10 logarithm of the permeability
figure('Position', [100 100, 900 260]); subplot(1,3,1)
plotCellData(G, lperm, 'edgecolor', 'none'); axis equal tight off
title('log10(K)','FontWeight','normal')
colormap(gca,flipud(pink))
% mrstColorbar(lperm,'south');

%% Filter 1
% We will make a density function based on the permeability as input to
% distmesh. To this end, we first construct a binary indicator (1 = high
% perm, 0 = low perm).
ind = lperm > 4;
% Next, we grow the indicator into the nearest neighbor a number of times
N  = getNeighbourship(G);
A = getConnectivityMatrix(N);
for j = 1:3
    ind = sum(A*ind,2) > 0;
end
ind = full(ind*1);
% Plot indicator
outlineCoarseGrid(G, ind, 'LineWidth',1);
title('Permeability + indicator','FontWeight','normal');

%% Filter 2
% To get a smooth transition from high- to low-density regions, we smoothen
% the indicator
S = 2*A + speye(G.cells.num);
S = S./sum(S,2);
for i = 1:20
    ind = S*ind;
end
% Plot indicator
subplot(1,3,2);cla
plotCellData(G, 1./(ind+.2), 'edgecolor', 'none'); axis equal tight off
title('Smoothed indicator','FontWeight','normal');
% mrstColorbar(1./(ind+.2),'south')
colormap(gca,flipud(parula))

%% Construct input functions to distmesh
% We use the inverse of the density (plus small number to avoid division by
% zero) as edge length function
x  = G.cells.centroids;
fh = scatteredInterpolant(x(:,1), x(:,2), 1./(ind + 0.2));
fh = @(x) fh(x(:,1), x(:,2));
% Grid input
L  = max(G.nodes.coords(:,1:2)); % Dimensions
n  = 80;                            % Approximate number of points
dx = max(L)/n;                      % Approximate cell diameter
% Distance function
rectangle = [0,0; L];
fd        = @(p,varargin) drectangle(p, 0, L(1), 0, L(2));
corners   = [0,0; 0,L(2); L(1),0; L(1),L(2)];
% Create voronoi sites with distmest
p = distmesh2d(fd, fh, dx, rectangle, corners, false);

%% Construct PEBI grid
% We construct a clipped PEBI grid using the voronoi sites
G1 = clippedPebi2D(p, [0, 0; L(1), 0; L(1), L(2); 0, L(2)]);
% Plot the grid
subplot(1,3,3)
plotCellData(G, lperm,'EdgeColor','none');
plotGrid(G1,'FaceColor','none','EdgeAlpha',.5); axis equal tight off
title('Adapted PEBI grid','FontWeight','normal');
colormap(gca,flipud(pink))
% mrstColorbar(lperm,'south')

%% Conforming grid
% Next, we construct a grid that conforms to the outlines of the
% hight-permeability channes. Again, we start with log10(perm_x) and apply
% filters
lperm = log10(modelRef.rock.perm(:,1));
lperm = lperm - min(lperm);

% Filter the permeability: replace by max value of all neighbors
nc = G.cells.num;
I   = sparse(N(:,1), N(:,2), max(lperm(N),[],2), nc, nc);
ind = full(max(max(I,[],2), max(I,[],1)'));

% Smoothen
for i = 1:3
    ind = S*ind;
end

% Find the channel outlines using contourc
xc = reshape(G.cells.centroids(:,1),G.cartDims);
yc = reshape(G.cells.centroids(:,2),G.cartDims);
c  = contourc(xc(:,1), yc(1,:)', reshape(ind, G.cartDims)',1);

% Extract the contour coordinates
[i, j] = deal(1);
permLines = cell(1,size(c,2));
d = max(L)*0.02;
while i < size(c,2)
    nPts = c(2,i);
    x = c(:, (1:nPts)+i)';
    x(x(:,1) < d | x(:,1) > L(1)-d | x(:,2) < d | x(:,2) > L(2)-d,:) = [];
    permLines{j} = x;
    i = i + nPts + 1;
    j = j+1;
end
permLines = permLines(1:j-1);

%% Plot the contours
figure('Position', [100 100, 1000 400]); subplot(1,2,1)
colormap(flipud(pink))
plotCellData(G, lperm, 'edgecolor', 'none');
plotLinePath(permLines,'-k','LineWidth',2);
axis equal tight off

%% Construct conforming grid
% First we find the outlines that are circular. pebiGrid expects to get
% lines, not circles, so we split these into two lines.
splitPermLines = {};
for i = 1:numel(permLines)
   if sum(abs(permLines{i}(1, :) - permLines{i}(end, :))) < 1e-10
      % Start and end point are equal -> circle.
      % Split line at midpoint
      num_points = size(permLines{i}, 1);
      l1 = permLines{i}(1 : ceil(num_points / 2), :);
      l2 = permLines{i}(ceil(num_points / 2) : end, :);
      splitPermLines{end + 1} = l1; %#ok<*SAGROW>
      splitPermLines{end + 1} = l2;
   else
      splitPermLines{end + 1} = permLines{i};
   end
end
% We construct a conforming grid using pebiGrid, with refinement around the
% outlines
dx = max(L)/10; % Grid cell size
G2 = pebiGrid2D(dx, L, ...
        'faceConstraints', permLines  , ... % Lines
        'interpolateFC'  , true       , ... % Interpolate faults
        'FCRefinement'   , true       , ... % Refine reservoir sites
        'FCFactor'       , 0.09       , ... % Relative fault cell size
        'FCEps'          , 0.07*max(L), ... % Size of ref transition region
        'linearize'      , true);
% Plot the grid
subplot(1,2,2)
plotCellData(G,lperm,'edgecolor','none')
plotGrid(G2,'facecolor','none','edgealpha',.5)
colormap(flipud(pink))
axis equal tight off

%% Sample rock properties
% We use sampleFromBox to sample petrophysical properties from the
% reference grid onto the two PEBI grids
G1 = computeGeometry(G1);
G2 = computeGeometry(G2);
% Sample permeability
dims  = G.cartDims(1:2);
perm1 = nan(G1.cells.num, G1.griddim);
perm2 = nan(G2.cells.num, G2.griddim);
for d = 1:G1.griddim
    perm1(:,d) = sampleFromBox(G1, reshape(rock.perm(:,d), G.cartDims));
    perm2(:,d) = sampleFromBox(G2, reshape(rock.perm(:,d), G.cartDims));
end
% Sample porosity
poro1 = sampleFromBox(G1, reshape(rock.poro, dims));
poro2 = sampleFromBox(G2, reshape(rock.poro, dims));
% Make rock structures
rock1 = makeRock(G1, perm1, poro1);
rock2 = makeRock(G2, perm2, poro2);

%%
figure('Position',[100 100 1120 420]);

subplot(1,3,1)
title('Density-based','FontWeight','normal') 
K = convertTo(rock1.perm(:,1),milli*darcy);
plotCellData(G1, log10(K), 'edgealpha', 0.2);
mrstColorbar(K,'south',true)
axis equal tight off

subplot(1,3,2)
title('Conforming','FontWeight','normal')
K = convertTo(rock2.perm(:,1),milli*darcy);
plotCellData(G2, log10(K), 'edgealpha', 0.2)
mrstColorbar(K,'south',true)
axis equal tight off

subplot(1,3,3)
title('Reference','FontWeight','normal')
K = convertTo(rock.perm(:,1),milli*darcy);
plotCellData(modelRef.G, log10(K), 'edgealpha', 0.2)
mrstColorbar(K,'south',true)
axis equal tight off
colormap(flipud(pink))

%% Simulate with conforming grid
% Finally, we simulate a pressure-drop problem with injection of water
% using the conforming PEBI grid
% Get top and bottom faces
faces  = boundaryFaces(G2);
tb     = abs(G2.faces.normals(faces,2)./G2.faces.areas(faces)) > 0.5;
bottom = tb & G2.faces.centroids(faces,2) < mean(G2.faces.centroids(faces,2));
top    = tb & G2.faces.centroids(faces,2) > mean(G2.faces.centroids(faces,2));
% add BC
bc = [];
bc = addBC(bc, faces(bottom), 'pressure', 100*barsa, 'sat', [1,0]);
bc = addBC(bc, faces(top), 'pressure', 1*barsa, 'sat', [1,0]);
% Schedule and initial state
time     = 2*year;
dt       = rampupTimesteps(time, 30*day);
schedule = simpleSchedule(dt, 'bc', bc);
state0   = initResSol(G2, 1*barsa, [0.2, 0.8]);
% Model
model = GenericBlackOilModel(G2, rock2, fluid, 'gas', false);
% Simulate
[~, states, rep] = simulateScheduleAD(state0, model, schedule);

%% Simulate with reference grid
% We compare the results with the reference grid
faces  = boundaryFaces(G);
tb     = abs(G.faces.normals(faces,2)./G.faces.areas(faces)) > 0.5;
bottom = tb & G.faces.centroids(faces,2) < mean(G.faces.centroids(faces,2));
top    = tb & G.faces.centroids(faces,2) > mean(G.faces.centroids(faces,2));
% Add BC
bcRef = [];
bcRef = addBC(bcRef, faces(bottom), 'pressure', 100*barsa, 'sat', [1,0]);
bcRef = addBC(bcRef, faces(top), 'pressure', 1*barsa, 'sat', [1,0]);
% Schedule and initial state
scheduleRef = simpleSchedule(dt, 'bc', bcRef);
state0Ref   = initResSol(G, 1*barsa, [0.2, 0.8]);
% Model
modelRef = GenericBlackOilModel(G, modelRef.rock, modelRef.fluid, 'gas', false);
% Simulate
[~, statesRef, repRef] = simulateScheduleAD(state0Ref, modelRef, scheduleRef);

%% Compare saturation profiles after time step 12
n = 12;
figure('Position', [100 100, 1000 400]);
bx1 = subplot(1,2,1);
bx1.Position(1)=bx1.Position(1)+.04;
plotCellData(G2,log10(rock2.perm(:,1)),'EdgeColor','none'); axis tight equal off
colormap(bx1,flipud(pink));
axlim = axis(bx1);
ax1 = axes('Position',bx1.Position);
plotCellData(G2,states{n}.s(:,1),states{n}.s(:,1)>.21,'EdgeColor','none');
axis equal tight off
colormap(ax1,flipud(winter));
axis(ax1,axlim), axis off

bx2 = subplot(1,2,2); 
plotCellData(G,log10(rock.perm(:,1)),'EdgeColor','none'); axis tight equal off
colormap(bx2,flipud(pink));
ax2 = axes('Position',bx2.Position);
plotCellData(G,statesRef{n}.s(:,1),statesRef{n}.s(:,1)>.21,'EdgeColor','none');
axis equal tight off
colormap(ax2,flipud(winter));
axis(ax2,axlim), axis off

% Alternative: for interactive visualization, uncomment the next two lines:
% spe10Fig('Conforming'), plotToolbar(G2, states)     ; axis equal tight
% spe10Fig('Reference') , plotToolbar(G, statesRef); axis equal tight

%% Compare simulation time
% We compare the simulation time for each timestep using the two grids
figure(), hold on
bar(1:numel(states), repRef.SimulationTime);
bar(1:numel(states), rep.SimulationTime, 0.4);
hold off; legend({'Reference grid', 'Conforming PEBI'});
