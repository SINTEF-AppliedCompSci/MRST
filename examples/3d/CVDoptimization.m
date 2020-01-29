%% This example shows how we can create a CVD grid. 
%We create a grid in the unit square. We also create a fault, and set these
%points as fixed points. 
%{
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2016-2020 Runar Lie Berge. See COPYRIGHT.TXT for details.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%}

%% set boundary to be the unit square
x = 1;
y = 1;
z = 1;
bdr   = [ 0, 0, 0;  ...
          x, 0, 0;  ...
          x, y, 0;  ...
          0, y, 0;  ...
          0, 0, z;  ...
          x, 0, z;  ...
          x, y, z;  ...
          0, y, z];

%% Set gridding parameters
fGs = x/5;
dt = x/5;
rho = @(p) fGs*(1+0*p(:,1));

%% Create Fault
swap1 = [0,1;1,0]==1;
swap2 = [1,0;0,1]==1;
rectangle1 = [min(bdr(:, [1,3])); max(bdr(:,[1,3]))];
fixedPts = [rectangle1; rectangle1(swap1)';rectangle1(swap2)'];
hd = @(p) rho(p)/fGs;
fd = @(p) drectangle(p, rectangle1(1),rectangle1(2), rectangle1(3),rectangle1(4));

[Pts,t] = distmesh2d(fd, hd, fGs, rectangle1, fixedPts, false);

fDt1.ConnectivityList = t;

fDt1.Points = [Pts(:,1), y/2*ones(size(Pts,1),1), Pts(:,2)];

F = surfaceSites3D({fDt1},{rho});

%% Create initial reservoir sites

xmax = max(bdr(:,1))-0.001; xmin = min(bdr(:,1))+0.001;
ymax = max(bdr(:,2))-0.001; ymin = min(bdr(:,2))+0.001;
zmax = max(bdr(:,3))-0.001; zmin = min(bdr(:,3))+0.001;


xr = xmin:dt:xmax;
yr = ymin:dt:ymax;
zr = zmin:dt:zmax;
[X,Y,Z] = ndgrid(xr,yr,zr);
rSites = [X(:), Y(:), Z(:)];
% Remove conflict sites
[rSites,removed] = surfaceSufCond3D(rSites,F.c.CC,F.c.R);

%% CVD optimization
bdrDT = delaunayTriangulation(bdr);
G = CPG3D(rSites,bdr,'tol',1e-3,'maxIt',20,'fixedPts',F.f.pts);
G = computeGeometry(G);

%% plot
% Grid view
figure() 
hold on
c = G.cells.centroids(:,3)>z*0.5 ...
  & G.cells.centroids(:,2)<y/2 ...
  & G.cells.centroids(:,1)>x/2; 
plotGrid(G,~c)
 xb = [0,x,x,0,0];
yb = [0,0,y,y,0];
zb = [0,0,0,0,0];

plot3(xb,yb,zb,'k')
plot3(xb,yb,zb+z,'k')
 for j = 1:numel(xb)
    plot3([xb(j),xb(j)],[yb(j),yb(j)],[0,z],'k')
 end

% make plot pretty
view(60,30);
light('Position',[-1 -1 -1],'Style','infinite')
light('Position',[-1 -1 2],'Style','local')
set(gca,'zdir','normal')
axis equal off
