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
fGs = x/8;
gs = x/8;
rho = @(p) fGs*ones(size(p,1),1);

%% Create Fault
swap1 = [0,1;1,0]==1;
swap2 = [1,0;0,1]==1;
rectangle1 = [min(bdr(:, [1,3])); max(bdr(:,[1,3]))];
fixedPts1 = [rectangle1; rectangle1(swap1)';rectangle1(swap2)'];
hd = @(p) rho(p)/fGs;
fd1 = @(p) drectangle(p, rectangle1(1),rectangle1(2), rectangle1(3),rectangle1(4));

[Pts1,t1] = distmesh2d(fd1, hd, fGs, rectangle1, fixedPts1);
fDt1.ConnectivityList = t1;

faultHeight1 = @(p) y/3*ones(size(p,1),1) + 0.3*p(:,1)+0.2*p(:,2);
fDt1.Points = [Pts1(:,1), faultHeight1(Pts1), Pts1(:,2)];

fDt = {fDt1};
F = createFaultGridPoints3D(fDt,{rho});

%% Create reservoir sites
n = round((1/gs)^3);
rSites = rand(n,3);

%% Enforce sufficient fault condition
nR = size(rSites,1);
rGs = zeros(nR,1);
rPri = zeros(nR,1);

[rSites,removed] = faultSufCond(rSites,F.c.CC,F.c.R);
rGs = rGs(~removed);
rPri = rPri(~removed);


%% Create Grid
pts = [F.f.pts;rSites];
gs  = [F.f.Gs;rGs];
pri = [F.f.pri;rPri];

bdrDT = delaunayTriangulation(bdr);
Gd = clippedPebi3D(pts,bdrDT);

%% Plot grid
plotGrid(Gd)
