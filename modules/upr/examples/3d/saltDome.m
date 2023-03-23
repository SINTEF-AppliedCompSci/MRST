% In this example we create a grid of a salt dome. We use the 3D distmesh
% routine to create a tetrahedrization of the salt dome, then we extract
% the boundary triangulation from this. 

%{
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2016-19 Runar Lie Berge. See COPYRIGHT.TXT for details.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%}

%% set boundary
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

fGs = x/10;
dt = x/10;
rho = @(p) fGs*ones(size(p,1),1);

%% Create fault
swap1 = [0,1;1,0]==1;
swap2 = [1,0;0,1]==1;

fixedPts = [];

%hd = @(p) 0.6*rho(p)/fGs;

x0 = [x/2, y/2,z/2];
h  = z/2.01;
r  = 0.3;
dsphr = @(p) dsphere(p,x0(1),x0(2),x0(3),r); % Sphere
dsyl  = @(p) dsylinder(p, x0,h,r);
ds    = @(p) dunion(dsphr(p), dsyl(p));

[p,t]=distmeshnd(ds,@huniform,fGs,[0,0,0;x,y,z],[]);
fDtVol = triangulation(t, p(:,1), p(:,2), p(:,3));
fDt.ConnectivityList = fDtVol.freeBoundary;
fDt.Points = p;
mid = 1/3*(fDt.Points(fDt.ConnectivityList(:,1),:) ...
         + fDt.Points(fDt.ConnectivityList(:,2),:) ...
         + fDt.Points(fDt.ConnectivityList(:,3),:));
fDt.ConnectivityList = fDt.ConnectivityList(mid(:,3)>x0(3)-h+1e-4,:);

F = surfaceSites3D({fDt},{rho});
clf
patch('vertices',fDt.Points, 'faces', fDt.ConnectivityList, ...
      'facealpha',0.3,'facecolor','r')
axis equal
%% Create reservoir sites

%rSites = rand(n,3);
xmax = max(bdr(:,1))-dt/2; xmin = min(bdr(:,1))+dt/2;
ymax = max(bdr(:,2))-dt/2; ymin = min(bdr(:,2))+dt/2;
zmax = max(bdr(:,3))-dt/2; zmin = min(bdr(:,3))+dt/2;


xg = xmin:dt:xmax;
yg = ymin:dt:ymax;
zg = zmin:dt:zmax;

[X,Y,Z] = ndgrid(xg,yg,zg);
rSites = [X(:), Y(:), Z(:)];

nR = size(rSites,1);
rGs = zeros(nR,1);
rPri = zeros(nR,1);

[rSites,removed] = surfaceSufCond3D(rSites,F.c.CC,F.c.R);
rGs = rGs(~removed);
rPri = rPri(~removed);

pts = [F.f.pts;rSites];
gs  = [F.f.Gs;rGs];
pri = [F.f.pri;rPri];

Gd = mirroredPebi3D(pts,bdr);

Gd = computeGeometry(Gd);
%% plot grid
figure(1);clf; hold on
plotGrid(Gd,'facecolor','none','edgealpha',0.4)
color = get(gca,'colororder');
patch('Vertices', fDt.Points, 'faces', fDt.ConnectivityList, ...
      'facecolor',color(2,:),'facealpha',1)

axis equal;
set(gca,'zdir','normal')
light('Position',[-.5,0,1])
view(-50,17)
axis([0,1,0,1,0,1])
axis equal tight off

figure(2); clf
cc = Gd.cells.centroids;

c1 = cc(:,3) < z/4.1;
c2 = cc(:,1) > x/2.1;
cf = ds(cc)<0;
plotGrid(Gd, (c1 | c2) & ~cf)
plotGrid(Gd, cf,'facecolor',color(2,:))
set(gca,'zdir','normal')
axis equal tight off
view(-50,17)
light('Position',[-.5,0,1])
