%% Section 4.3:  PEBI conformity
% We show an example with three intersection fracture planes to demonstrate
% how to build a hierarchy of conforming PEBI grid, starting from the 1D
% intersections of each pair of planes, each of which is tessellated by a
% 2D grid, and ending up in the 3D volumetric grid.

%% Define the fracture planes as polygons and compute intersection
poly1 = ellipticSurface3D([3,3, 3], 1.5, 1.5, 0, 0, pi/2);
poly2 = [2,2,3.3; 5,2,3.3; 5,4,3.3; 2,4, 3.3];
poly3 = [3,1,1; 3,5,1; 3,5,5; 3,1, 5];
% Collect fractures
fracs = {poly1, poly2, poly3};

% Plot the fractures
fig1 = figure();
col  = get(gca, 'colororder');
grey = [0.9, 0.9, 0.9];
hf   = zeros(numel(fracs),1);
hold on
for i = 1:numel(fracs)
    f = fracs{i};
    hf(i) = patch(f(:,1), f(:,2), f(:,3), ones(size(f,1),1), ...
        'LineWidth',1,'FaceAlpha',.2,'FaceColor',col(i,:)); 
end
view(3)

%% Calculate intersections and tessellate them (1D grids)
% Calculate intersections
intersections = surfaceIntersections3D(fracs);

% Tessellate the intersections
ds      = 0.25;            % mesh size  
gamma   = ds./[1 4 8];     % fracture size distance
grids1D = lineGrid3D(intersections, gamma(1));

figure(fig1);
hp = zeros(numel(intersections),1);
for i = 1:numel(intersections)
   int = intersections{i}{1}; 
   pts = grids1D{i}{1};
   plot3(int(:,1), int(:,2), int(:,3))
   hp(i) = plot3(pts(:,1), pts(:,2), pts(:,3),'.', 'markersize',15);
end

%%
% Plot the 1D grids that are witin the elliptic fault. Before plotting, we
% translate and rotate 
fig2   = figure(); hold on
center = mean(fracs{1}, 1);
R = rotationMatrixFromPlane(bsxfun(@minus, fracs{1}, center));
for k = 1:numel(grids1D)
    sites = bsxfun(@minus, grids1D{k}{1}, center)*R';
    nodes = bsxfun(@minus, grids1D{k}{2}, center)*R';
    if isempty(nodes) || any(abs(nodes(:,3))>1e-6)
        continue
    end
    plot(nodes(:, 1), nodes(:, 2), 'color', col(k, :), 'linewidth', 3)
    plot(sites(:, 1), sites(:, 2), '.m','markersize', 20)
end
set(gca, 'zdir', 'reverse')
axis tight equal off

%% Tessellate the fracture surfaces (2D grids)
grids2D = surfaceGrid3D(fracs, grids1D, intersections, ds, gamma(2));

figure(fig1);
set(hf,'FaceAlpha',1);
set(hp,'MarkerSize',20,'Color','k');
for i = 1:numel(grids2D)
    pts = grids2D{i}.cells.sites;
    plot3(pts(:,1), pts(:,2), pts(:,3),'.', 'markersize',5)
end

%%
% Show the tessellated elliptic fault, as in Figure 8
figure(fig2)
G = grids2D{1};
G.nodes.coords = bsxfun(@minus, G.nodes.coords, center) * R';
G.nodes.coords(:, 3) = 0.01;
G.cells.sites = bsxfun(@minus, G.cells.sites, center);
G.cells.sites = G.cells.sites * R';
plotGrid(G, 'faceColor', grey);
i = G.cells.resSite;
plot(G.cells.sites( i,1), G.cells.sites( i,2), 'k.', 'markersize',10)
plot(G.cells.sites(~i,1), G.cells.sites(~i,2), '.g', 'markersize',10);

%%
% Plot outline of the 2D and 1D grids in 3D as in the left plot of Figure 8
figure
col = get(gca, 'colororder');
for k =1:numel(grids2D)
    plotGrid(grids2D{k}, 'facecolor', grey, 'edgealpha', 0)
    bnd_faces = any(grids2D{k}.faces.neighbors == 0, 2);
    plotFaces(grids2D{k}, bnd_faces)
end
for k = 1:numel(grids1D)
    g1 = grids1D{k}{2};
    plot3t(g1(:, 1), g1(:,2), g1(:,3), 0.05, col(k, :))%, 'linewidth', 6)
end
set(gca, 'zdir', 'normal'), view(160, 20)
axis equal off tight,  light('position', [10,5,10])

%% Define the 3D volumetric grid
G3 = volumeGrid3D([6, 6, 6], fracs, grids2D, ds, gamma(3));
G3 = computeGeometry(G3);
G3 = mrstGridWithFullMappings(G3);

%% Plot the 1D and 2D grids
% These grids are defined in different coordinate systems, so we must
% rotate and translate them to get the plots shown in the upper part of
% Figure 9
figure('position',[100 100 980 420])
for i = 1:numel(grids2D)
    G = grids2D{i};
    center = mean(fracs{i}, 1);
    R  = rotationMatrixFromPlane(bsxfun(@minus, fracs{i}, center));
    G.nodes.coords = bsxfun(@minus, G.nodes.coords, center) * R';
    G.nodes.coords(:, 3) = 0.01;
    subplot(1,3,i)
    set(gca, 'zdir', 'normal')
    hold on
    plotGrid(G, 'faceColor', grey);%col(7 - 2*i,:))
    for k = 1:numel(grids1D)
        g1 = grids1D{k}{2};
        if isempty(g1)
            continue
        end
        nodes = bsxfun(@minus, g1, center);
        nodes = nodes * R';
        if abs(nodes(1, 3)) > 1e-6
            continue
        end
        plot(nodes(:, 1), nodes(:, 2), 'color', col(k, :), 'linewidth', 3)
    end
    axis tight equal off
end

%% Plot the 3D grid
% To visualize the inside, we split it in two along the plane defined by
% the hydraulic fracture (the elliptic polygonal) to get the plot shown in
% the lower part of Figure 9

% One half of the grid
figure('position',[100 100 980 420]);
subplot(1,2,1)

% Plot the 1D grids for the intersections
n  = normalFromPoints(poly1);
for k = 1:numel(grids1D)
    g1 = grids1D{k}{2};
    if isempty(g1) || ~all(bsxfun(@minus, g1, poly1(1,:))* n' > -1e-3)
        continue
    end
    plot3t(g1(:, 1), g1(:,2), g1(:,3), 0.05, col(k, :))
end

% Plot the half of the grid
c  = true(G3.cells.num,1);                   
c  = c & (bsxfun(@minus, G3.cells.centroids, poly1(1,:))* n' > 1e-5);
plotGrid(G3, c)

% Plot the gridded hydraulic fracture in gray
ffaces = (abs(G3.faces.centroids(:,2) - 3) < 1e-5) ...
    & (abs(G3.faces.normals(:, 2))./G3.faces.areas > 1-1e-6);
plotFaces(G3, ffaces, 'facecolor', grey);

set(gca, 'zdir', 'normal')
view(130, 10)
axis equal off tight
light('position', [3,10,10])

% The second half of the grid
subplot(1,2,2), cla
plotGrid(G3,~c);
plotFaces(G3, ffaces, 'facecolor', grey);

n  = -normalFromPoints(poly1);
for k = 1:numel(grids1D)
    g1 = grids1D{k}{2};
    if isempty(g1) || ~all(bsxfun(@minus, g1, poly1(1,:))* n' > -1e-3)
        continue
    end
    plot3t(g1(:, 1), g1(:,2), g1(:,3), 0.05, col(k, :))
end

set(gca, 'zdir', 'normal')
view(45, 10)
axis equal off tight
light('position', [3,-10,-10])