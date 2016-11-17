%% Grid-orientation Errors for Five-Spot
% Compare solutions on a quarter five-spot and a rotated five-spot
mrstModule add incomp


%% Fluid model
T = 5*year; nstep=1;
gravity reset off
fluid = initSimpleFluid('mu' , [   1,     1] .* centi*poise     , ...
                        'rho', [1000,  850] .* kilogram/meter^3, ...
                        'n'  , [   2,    2]);

%% Parallel grid
Gp = computeGeometry(cartGrid([20,20,1],[1000 1000 20]));
prock = makeRock(Gp, 450*milli*darcy, 0.2);
rate  = .6*sum(poreVolume(Gp,prock))/T;
Wp = addWell([],Gp, prock, 1, 'Type', 'rate', ...
    'Val', rate, 'name', 'I', 'radius', .1, 'Comp_i', [1 0]);
Wp = addWell(Wp,Gp, prock, Gp.cells.num, 'Type', 'bhp', ...
    'Val', -rate, 'name', 'P', 'radius', .1, 'Comp_i', [0 1]);

xp  = initState(Gp,Wp,100*barsa, [0 1]);
hTp = computeTrans(Gp, prock);
for i=1:nstep
    xp  = incompTPFA(xp, Gp, hTp, fluid, 'wells', Wp);
    xp  = explicitTransport(xp, Gp, T/nstep, prock, fluid, 'wells', Wp);
end
plotCellData(Gp,xp.s(:,1));
plotWell(Gp, Wp);
%x = Gp.cells.centroids(vertcat(Wp.cells),:);
%hold on; plot(x(:,1),x(:,2),'ow'); hold off;

%% Diagonal grid
theta = pi/4;
R = [cos(theta) -sin(theta); sin(theta) cos(theta)];
Gd = cartGrid([round(Gp.cartDims(1:2).*sqrt(2)) 1],[1000 1000 20]);
Gd.nodes.coords(:,1:2) = sqrt(2)*(R*(Gd.nodes.coords(:,1:2)'))';
Gd = computeGeometry(Gd);
drock = makeRock(Gd, 450*milli*darcy, 0.2);
n = Gd.cartDims(1);
Wd = addWell([], Gd, drock, 1, 'Type', 'rate', ...
    'Val', rate, 'name', 'I', 'radius', .1, 'Comp_i', [1 0]);
Wd = addWell(Wd, Gd, drock, n*n, 'Type', 'rate', ...
    'Val', rate, 'name', 'I', 'radius', .1, 'Comp_i', [1 0]);
Wd = addWell(Wd, Gd, drock, n, 'Type', 'rate', ...
    'Val', -rate, 'name', 'P', 'radius', .1, 'Comp_i', [0 1]);
Wd = addWell(Wd, Gd, drock, n*(n-1)+1, 'Type', 'rate', ...
    'Val', -rate, 'name', 'P', 'radius', .1, 'Comp_i', [0 1]);

xd  = initState(Gd, Wd, 100*barsa, [0 1]);
hTd = computeTrans(Gd, drock);
for i=1:nstep
    xd  = incompTPFA(xd, Gd, hTd, fluid, 'wells', Wd);
    xd  = explicitTransport(xd, Gd, T/nstep, drock, fluid, 'wells', Wd);
end
plotCellData(Gd,xd.s(:,1));
plotWell(Gd,Wd);
axis equal
%x = Gd.cells.centroids(vertcat(Wd.cells),:);
%hold on; plot(x(:,1),x(:,2),'*k'); hold off;

%% Remap the diagonal grid for plotting
rp = xp.s(:,1);
rd = xd.s(:,1);
[i,j] = gridLogicalIndices(Gd);
n = max(i);
m = round(n/2);
k = mod(n,2);
map = [m+1:n 1:m]';
ind = sub2ind(Gd.cartDims,i,map(j));
rmap = rd(ind);
rmap = reshape(rmap,Gd.cartDims);
rmap(1:end,1:m-k) = rmap(end:-1:1,1:m-k);
clf
subplot(1,2,1);
plotCellData(Gd,rd); axis equal tight
subplot(1,2,2);
mx = max(Gp.nodes.coords(:,1:2));
Gd.nodes.coords(:,1:2) = bsxfun(@plus, Gd.nodes.coords(:,1:2),mx.*[.5 -.5]);
Gd = computeGeometry(Gd);
plotCellData(Gd,rmap(:)); axis equal tight

%%
clf
cval = linspace(0,1,11); cval=.5*cval(1:end-1)+.5*cval(2:end);
contour(reshape(Gp.cells.centroids(:,1), Gp.cartDims),...
    reshape(Gp.cells.centroids(:,2), Gp.cartDims), ...
    reshape(rp,Gp.cartDims), ...
    cval, '--','LineWidth',3);
hold on
contour(reshape(Gd.cells.centroids(:,1), Gd.cartDims),...
    reshape(Gd.cells.centroids(:,2), Gd.cartDims), ...
    reshape(rmap(:),Gd.cartDims)+5e-3*rand(Gd.cartDims), ...
    cval, '-','LineWidth',.5);
hold off
axis equal tight; axis([0 mx(1) 0 mx(2)]);