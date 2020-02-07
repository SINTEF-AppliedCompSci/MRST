%% Grid adapting to a multilateral well
% We consider a multilateral well with four branches. Here, the branches
% are assumed to be in a single vertical plane to simplify visualization of
% the resulting grid, but the procedure we present is more general and can
% easily be changed to handle horizontal well paths.

%% Set gridding parameters
xmax = 1.001; xmin = -0.001;
ymax = 1.001; ymin = -0.001;
zmax = 1.001; zmin = -0.001;
bdr   = [ xmin, ymin, zmin;  ...
          xmax, ymin, zmin;  ...
          xmax, ymax, zmin;  ...
          xmin, ymax, zmin;  ...
          xmin, ymin, zmax;  ...
          xmax, ymin, zmax;  ...
          xmax, ymax, zmax;  ...
          ymin, ymax, zmax];

%% Create well sites
wGc = @(p) 1/25 / 2 * ones(size(p, 1), 1);
line1 = [0.5, 0.5, 1; 0.5, 0.5, 3 / 4];
line2 = [0.5, 0.5, 3 / 4; 1 / 5, 0.5, 0.5; 1 / 6, 0.5, 1 / 6];
line3 = [0.5, 0.5, 3 / 4; 4 / 5, 0.5, 0.5; 5 / 6, 0.5, 1 / 6];
line4 = [0.5, 0.5, 3 / 4; 0.5, 0.5, 0.5];
line5 = [0.5, 0.5, 0.5; 3 / 5, 0.5, 1 / 3; 2 / 3, 0.5, 0];
line6 = [0.5, 0.5, 0.5; 2 / 5, 0.5, 1 / 3; 1 / 3, 0.5, 0];
wl = {line1, line2, line3, line4, line5, line6};

figure, hold all
for i=1:numel(wl)
    plot3(wl{i}(:,1),wl{i}(:,2),wl{i}(:,3),'LineWidth',2);
end
view(3), axis tight off
zoom(1.5), set(gca,'Clipping',false)
%%
Wc = lineSites3D(wl, repmat({wGc}, 1, numel(wl)));
plot3(Wc.pts(:,1),Wc.pts(:,2),Wc.pts(:,3),'.','MarkerSize',10);


%% Create reservoir sites
dt = 1 / 13; % Mesh size
xr = xmin + dt / 2 : dt : xmax - dt / 2;
yr = ymin + dt / 2 : dt : ymax - dt / 2;
zr = zmin + dt / 2 : dt : zmax - dt / 2;
[X, Y, Z] = ndgrid(xr, yr, zr);
bgSites = [X(:), Y(:), Z(:)];

varArg = {'level', 1, 'maxLev', 2, 'distTol', 1.5 * dt};
res = {};
for i = 1:size(bgSites, 1)
    res = [res; mlqt(bgSites(i,:), Wc.pts , [dt, dt, dt], varArg{:})];
end
bgSites = vertcat(res{:, 1});

% Remove any sites innside the spheres
bgSites = lineSufCond3D(bgSites, Wc);

% Plot result
plot3(bgSites(:,1),bgSites(:,2),bgSites(:,3),'.');
plotGrid(cartGrid([1 1 1],[xmax,ymax,zmax]),'FaceColor','none');
set(gca,'zdir','normal')
axis tight

%% Create the 3D grid
sites  = [Wc.pts; bgSites];
G      = mirroredPebi3D(sites, bdr);

%% Plots of the grid
% Plot the the perforated cells only
G = computeGeometry(G);
color = get(gca, 'colororder');
figure()
plotGrid(G, 1:size(Wc.pts, 1), 'facecolor', color(1, :))
axis equal off tight
set(gca,'zdir','normal')
view(0,0)

%%
% Plot the grid split open in two
figure(2); clf
cl = sites(:, 2) < 0.5 + 1e-5;
cr = sites(:, 2) >= 0.5 + 1e-5;
Gp = G;
plotGrid(Gp, cl)
plotGrid(Gp, 1:size(Wc.pts, 1), 'facecolor', color(1, :))
Gp = rodriguesRotation(Gp, -pi / 2, [0, 0, 1], [1, 0.5, 0]);
Gp.nodes.coords(:, 2) = Gp.nodes.coords(:, 2) + 0.2;
plotGrid(Gp, cr)
view(-130, 25)
light('Position',[0.5 1.5 1.5],'Style','local')
light('Position',[0.5 1.5 1.5],'Style','infinite')
axis equal tight off
set(gca,'zdir','normal')

%%
% Plot half of the grid, seen from another direction
figure(3); clf
ch = sites(:, 1) < 0.5 + 1e-5;
Gp = G;
plotGrid(Gp, ch)
plotGrid(Gp, 1:size(Wc.pts, 1), 'facecolor', color(1, :))
view(140, 30)
light('Position',[1 1.5 1.5],'Style','local')
light('Position',[1 1.5 1.5],'Style','infinite')
axis equal tight off
set(gca,'zdir','normal')