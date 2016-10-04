%% Set boundary to the unit cube
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
wGs = x/8;
gs = x/8;
rho = @(p) wGs*ones(size(p,1),1);

%% Create well
wellLine = [0.5,0.5,0.01;0.8,0.5,0.5;0.8,0.5,0.99];
W        = createWellGridPoints3D({wellLine},rho);
%% Create reservoir sites
n = round((1/gs)^3);
rSites = rand(n,3);

%% Enforce sufficient well condition
[rSites,removed] = wellSufCond3D(rSites, W);


%% Create Grid
pts = [W.pts;rSites];

bdrDT = delaunayTriangulation(bdr);
Gd = clippedPebi3D(pts,bdrDT);
Gd = computeGeometry(Gd);
%% Plot grid
plotGrid(Gd)
plotGrid(Gd,1:size(wellPts,1),'facecolor','r')
axis equal
figure()
c = Gd.cells.centroids(:,2)<0.5;
plotGrid(Gd,c)
plotGrid(Gd,1:size(W.pts,1),'facecolor','r')
axis equal