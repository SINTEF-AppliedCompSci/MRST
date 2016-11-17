%% Grid-orientation effects for the five-spot
% With a two-point type discretization (and with many other discretizations
% as well), evolving displacement profiles will preferrentially move along
% the axial directions, i.e., in the direction of the normals to the cell
% faces. To illustrate this, we contrast approximate solutions for the
% repeated five-spot well pattern computed using the original quarter
% five-spot setup and with a rotated setup with injectors in the SW and NE
% corners and producers in the SE and NW corners. The two setups will have
% different preferrential flow directions and hence generally give
% different solutions.

mrstModule add incomp

%% Common parts of the two simulation setups
T = 5*year; 
% nstep=1;
% cartDim = [16 16 1];
pvi = 0.6;
gravity reset off
fluid = initSimpleFluid('mu' , [   1,     1] .* centi*poise     , ...
                        'rho', [1000,  850] .* kilogram/meter^3, ...
                        'n'  , [   2,    2]);

%% Original quater-five spot
% This is a standard quarter-five spot with injector and producer in
% diagonally opposit corners. Here the main flow direction will be along
% the diagonal of the grid cells.
domain = [1000 1000 20];
G      = computeGeometry(cartGrid(cartDim,domain));
rock   = makeRock(G, 450*milli*darcy, 0.2);
rate   = pvi*sum(poreVolume(G,rock))/T;
W = addWell([],G, rock, 1, 'Type', 'rate', ...
    'Val', rate, 'name', 'I', 'radius', .1, 'Comp_i', [1 0]);
W = addWell(W,G, rock, G.cells.num, 'Type', 'rate', ...
    'Val', -rate, 'name', 'P', 'radius', .1, 'Comp_i', [0 1]);

x  = initState(G,W,100*barsa, [0 1]);
hT = computeTrans(G, rock);
for i=1:nstep
    x  = incompTPFA(x, G, hT, fluid, 'wells', W);
    x  = explicitTransport(x, G, T/nstep, rock, fluid, 'wells', W);
end
%{
clf
plotCellData(G,x.s(:,1));
plotWell(G, W);
%}

% Set plotting function for later
N = 10;
cval = linspace(0,1,N+1); cval=.5*cval(1:end-1)+.5*cval(2:end);
plotData = @(G,x) contour(reshape(G.cells.centroids(:,1),G.cartDims), ...
                          reshape(G.cells.centroids(:,2),G.cartDims), ...
                          reshape(x,G.cartDims), cval, '-k', 'LineWidth',1.5);
colormap(flipud(.5*jet(N)+.5*ones(N,3)));

%% Rotated quater-five spot
% This is a grid that is rotated 45 degrees so that we simulate a case with
% four wells with injectors in the SW and NE corners and producers in the
% SE and NW corners. Here the main flow direction will align with the axial
% directions of the grid.
theta = pi/4;
R = [cos(theta) -sin(theta); sin(theta) cos(theta)];
Gr = cartGrid([round(G.cartDims(1:2).*sqrt(2)) 1], domain);
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
%{
clf
plotCellData(Gr,xr.s(:,1),'EdgeColor',[.6 .6 .6]);
plotWell(Gr,Wr);
i = G.cells.centroids(:,1)>=G.cells.centroids(:,2);
plotCellData(G,x.s(:,1),i,'EdgeColor','k');
plotWell(G, W);
axis tight equal;
c = jet; colormap([1 1 1; c(2:end,:)]);
%}

%% Compare the two solutions
remapAndPlot(Gr,xr.s(:,1),domain([1 2]),cval);
hold on; plotData(G,x.s(:,1)); hold off;
axis equal tight; axis([0 domain(1) 0 domain(2)]);
%legend('Original','Rotated',4);