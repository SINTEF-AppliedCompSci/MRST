%% A 3D domain with 2 curved faults. 
% We create two curved faults by triangulating them by distmesh.
% We then use the createFaultGridPoint3d function to create the grid points
% of the faults, before we crate a cartesian background grid.
%{
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2016-2020 Runar Lie Berge. See COPYRIGHT.TXT for details.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%}
%% set boundary
tx = 1;
ty = 1;
tz = .4;
bdr   = [ 0, 0, 0;  ...
          tx, 0, 0;  ...
          tx, ty, 0;  ...
          0, ty, 0;  ...
          0, 0, tz;  ...
          tx, 0, tz;  ...
          tx, ty, tz;  ...
          0, ty, tz];

%% Set gridding parameters
fGs = tx/20;
rho = @(p) fGs*(1+0*p(:,1));

%% Triangulate surfaces
swap1 = [0,1;1,0]==1;
swap2 = [1,0;0,1]==1;
rectangle1 = [min(bdr(:, [1,3])); max(bdr(:,[1,3]))];
fixedPts = [rectangle1; rectangle1(swap1)'; rectangle1(swap2)'];
hd = @(p) rho(p)/fGs;
fd = @(p) drectangle(p, rectangle1(1),rectangle1(2), rectangle1(3), rectangle1(4));

[Pts,t] = distmesh2d(fd, hd, fGs, rectangle1, fixedPts, false);

t1.ConnectivityList = t;
t2.ConnectivityList = t;

faultHeight1z = @(p)   ty/6*ones(size(p,1),1) + 0.4*p(:,2);
faultHeight2z = @(p) 5*ty/6*ones(size(p,1),1) - 0.4*p(:,2);

t1.Points = [Pts(:,1), faultHeight1z(Pts), Pts(:,2)];
t2.Points = [Pts(:,1), faultHeight2z(Pts), Pts(:,2)];

faultHeight1x = @(p) p(:,2) + 0.2*(p(:,1)-0.75).^2 + 0.1*sin(2*pi/tx*p(:,1));
faultHeight2x = @(p) p(:,2) + 0.2*(p(:,1)-0.75).^2;
t1.Points(:,2) = faultHeight1x(t1.Points);
t2.Points(:,2) = faultHeight2x(t2.Points);

figure
color=get(gca,'colororder');
hold on
trisurf(t1.ConnectivityList,t1.Points(:,1),t1.Points(:,2),t1.Points(:,3),...
    'FaceColor',color(1,:),'EdgeAlpha',.1);
trisurf(t2.ConnectivityList,t2.Points(:,1),t2.Points(:,2),t2.Points(:,3),...
    'FaceColor',color(2,:),'EdgeAlpha',.1);
view(3);

%% Place sites at sphere intersections
% t1 and t2 are triangulations of the two surfaces
R = @(p) 1/20 * ones(size(p, 1), 1); % Radius of spheres
F = surfaceSites3D({t1, t2}, {R, R});

plot3(F.f.pts(:,1),F.f.pts(:,2),F.f.pts(:,3),'.',...
    'MarkerSize',14,'Color',color(3,:));
plotGrid(cartGrid([1 1 1],[tx ty tz]),'FaceColor','none');
set(gca,'Projection','perspective')
view(-100,50); camlight left, camlight headlight

bdr = 1.01 * bdr;

%% Create reservoir sites
dt = 1 / 20; % Mesh size
xmax = max(bdr(:,1)) - dt / 2; xmin = min(bdr(:,1)) + dt / 2;
ymax = max(bdr(:,2)) - dt / 2; ymin = min(bdr(:,2)) + dt / 2;
zmax = max(bdr(:,3)) - dt / 2; zmin = min(bdr(:,3)) + dt / 2;

xr = xmin:dt:xmax; yr = ymin:dt:ymax; zr = zmin:dt:zmax;
[X,Y,Z] = ndgrid(xr, yr, zr);
rSites = [X(:), Y(:), Z(:)];
% Remove any sites innside the spheres
rSites = surfaceSufCond3D(rSites, F.c.CC, F.c.R);

sites = [F.f.pts; rSites];
G = mirroredPebi3D(sites, bdr);
plot3(X(:),Y(:),Z(:),'.','MarkerSize',6,'Color',color(4,:));
%% plot
G = computeGeometry(G);
color = get(gca,'ColorOrder');

% Grid view
figure(); hold on
cl = G.cells.centroids(:,2)>faultHeight2z(G.cells.centroids(:,[1,3])) + ...
                           faultHeight2x(G.cells.centroids(:,[1,3])) - ...
                           G.cells.centroids(:,3);
cr = G.cells.centroids(:,2)<faultHeight1z(G.cells.centroids(:,[1,3])) + ...
                           faultHeight1x(G.cells.centroids(:,[1,3])) - ...
                           G.cells.centroids(:,3);
                  
cm = ~cl & ~cr;

plotGrid(G,cl)

% Shift grid
dy = 0.3;
Gs = G;
Gs.nodes.coords(:, 2) = Gs.nodes.coords(:, 2) - dy;
plotGrid(Gs, cm)
Gs.nodes.coords(:, 2) = Gs.nodes.coords(:, 2) - dy;
plotGrid(Gs, cr)

tx = t1.Points(:, 1);
ty = t1.Points(:, 2) - 1.5 * dy;
tz = t1.Points(:, 3);
trisurf(t1.ConnectivityList, tx, ty, tz,'facecolor',color(1, :))

tx = t2.Points(:, 1);
ty = t2.Points(:, 2) - .5 * dy;
tz = t2.Points(:, 3);
trisurf(t2.ConnectivityList, tx, ty, tz,'facecolor',color(2, :))

% make plot pretty
view(120,30);
light('Position',[5 10 -1],'Style','infinite')
light('Position',[-1 -1 2],'Style','local')
set(gca,'zdir','normal')
axis equal off
