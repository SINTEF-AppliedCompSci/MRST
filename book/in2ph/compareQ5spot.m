%% Constants 
T = 5*year; nstep=128;
cartDim = [64 64 1];

%% Two-phase model
gravity reset off
fluid = initSimpleFluid('mu' , [   1,     1] .* centi*poise     , ...
                        'rho', [1000,  850] .* kilogram/meter^3, ...
                        'n'  , [   2,    2]);

%% Original quater-five spot
% This is a standard quarter-five spot with injector and producer in
% diagonally opposit corners. Here the main flow direction will be
domain = [1000 1000 20];
G      = computeGeometry(cartGrid(cartDim,domain));
rock   = makeRock(G, 450*milli*darcy, 0.2);
rate   = .6*sum(poreVolume(G,rock))/T;
W = addWell([],G, rock, 1, 'Type', 'rate', ...
    'Val', rate, 'name', 'I', 'radius', .1, 'Comp_i', [1 0]);
W = addWell(W,G, rock, G.cells.num, 'Type', 'bhp', ...
    'Val', -rate, 'name', 'P', 'radius', .1, 'Comp_i', [0 1]);

x  = initState(G,W,100*barsa, [0 1]);
hT = computeTrans(G, rock);
for i=1:nstep
    x  = incompTPFA(x, G, hT, fluid, 'wells', W);
    x  = explicitTransport(x, G, T/nstep, rock, fluid, 'wells', W);
end
clf
plotCellData(G,x.s(:,1));
plotWell(G, W);

%% Rotated quater-five spot
% This is a grid that is rotated 45 degrees so that we simulate a case with
% four wells with injectors in the SW and NE corners and producers in the
% SE and NW corners.
theta = pi/4;
R = [cos(theta) -sin(theta); sin(theta) cos(theta)];
Gr = cartGrid([round(G.cartDims(1:2).*sqrt(2)) 1],[1000 1000 20]);
Gr.nodes.coords(:,1:2) = sqrt(2)*(R*(Gr.nodes.coords(:,1:2)'))';
Gr = computeGeometry(Gr);
rockr = makeRock(Gr, 450*milli*darcy, 0.2);
n = Gr.cartDims(1);
Wr = addWell([], Gr, rockr, 1, 'Type', 'rate', ...
    'Val', rate, 'name', 'I', 'radius', .1, 'Comp_i', [1 0]);
Wr = addWell(Wr, Gr, rockr, n*n, 'Type', 'rate', ...
    'Val', rate, 'name', 'I', 'radius', .1, 'Comp_i', [1 0]);
Wr = addWell(Wr, Gr, rockr, n, 'Type', 'rate', ...
    'Val', -rate, 'name', 'P', 'radius', .1, 'Comp_i', [0 1]);
Wr = addWell(Wr, Gr, rockr, n*(n-1)+1, 'Type', 'rate', ...
    'Val', -rate, 'name', 'P', 'radius', .1, 'Comp_i', [0 1]);

xr  = initState(Gr, Wr, 100*barsa, [0 1]);
hTr = computeTrans(Gr, rockr);
for i=1:nstep
    xr  = incompTPFA(xr, Gr, hTr, fluid, 'wells', Wr);
    xr  = explicitTransport(xr, Gr, T/nstep, rockr, fluid, 'wells', Wr);
end
clf
plotCellData(Gr,xr.s(:,1),'EdgeColor',[.6 .6 .6]);
plotWell(Gr,Wr);
i = G.cells.centroids(:,1)>=G.cells.centroids(:,2);
plotCellData(G,x.s(:,1),i,'EdgeColor','k');
plotWell(G, W);
axis tight equal;
c = jet; colormap([1 1 1; c(2:end,:)]);

%% Remap the diagonal grid for plotting
rp = x.s(:,1);
rd = xr.s(:,1);
[i,j] = gridLogicalIndices(Gr);
n = max(i);
m = round(n/2);
k = mod(n,2);
map = [m+1:n 1:m]';
ind = sub2ind(Gr.cartDims,i,map(j));
rmap = rd(ind);
rmap = reshape(rmap,Gr.cartDims);
rmap(1:end,1:m-k) = rmap(end:-1:1,1:m-k);
clf
subplot(1,2,1);
plotCellData(Gr,rd); axis equal tight
hold on; plot([0 0 1000 1000 0],[0 1000 1000 0 0],'k','LineWidth',2); hold off
subplot(1,2,2);
mx = max(G.nodes.coords(:,1:2));
Gr.nodes.coords(:,1:2) = bsxfun(@plus, Gr.nodes.coords(:,1:2),mx.*[.5 -.5]);
Gr = computeGeometry(Gr);
plotCellData(Gr,rmap(:)); axis equal tight
hold on; plot([0 0 mx([1 1]) 0],[0 mx([1 1]) 0 0],'w','LineWidth',2); hold off

%%
clf
cval = linspace(0,1,11); cval=.5*cval(1:end-1)+.5*cval(2:end);
contour(reshape(G.cells.centroids(:,1), G.cartDims),...
    reshape(G.cells.centroids(:,2), G.cartDims), ...
    reshape(rp,G.cartDims), ...
    cval, '--','LineWidth',3);
hold on
contour(reshape(Gr.cells.centroids(:,1), Gr.cartDims),...
    reshape(Gr.cells.centroids(:,2), Gr.cartDims), ...
    reshape(rmap(:),Gr.cartDims)+5e-3*rand(Gr.cartDims), ...
    cval, '-','LineWidth',.5);
hold off
axis equal tight; axis([0 mx(1) 0 mx(2)]);
legend('Original','Rotated',4);