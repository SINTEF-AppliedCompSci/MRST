function runQ5DiagParal(cartDim, nstep, fluid, pvi)
%Solve and plot quarter five-spot problem on diagonal and parallel grid
%
% SYNOPSIS:
%   runQ5DiagParal(nstep, cartDim, fluid, pvi)
%
% DESCRIPTION:
%   This function runs two quarter five-spot problems and plots them in the
%   same axes. The first is the original quarter-five spot with injector
%   and producer in diagonally opposite corners. Here, the main flow
%   direction will be along the diagonal of the grid cells. This solution
%   will be plotted using contour lines. The second has a grid that is
%   rotated 45 degrees so that we simulate a case with four wells with
%   injectors in the SW and NE corners and producers in the SE and NW
%   corners. Here the main flow direction will align with the axial
%   directions of the grid. This is plotted as filled contours.
%
% REQUIRED PARAMETERS:
%   cartDim - Dimension of Cartesian grid used for the original Q5 setup
%
%   nstep   - Number of time steps in sequential solution
%
%   fluid   - Fluid object as defined by function 'initSimpleFluid'
%
%   pvi     - Dimensionless end time for simulation 

gravity reset off


% -------------------------------------------------------------------------
% Original quarter-five spot
% This is a standard 
domain = [1000 1000 20];
G      = computeGeometry(cartGrid(cartDim,domain));
rock   = makeRock(G, 450*milli*darcy, 0.2);
rate   = pvi*sum(poreVolume(G,rock));
W = addWell([],G, rock, 1, 'Type', 'rate', ...
    'Val', rate, 'name', 'I', 'radius', .1, 'Comp_i', [1 0]);
W = addWell(W,G, rock, G.cells.num, 'Type', 'rate', ...
    'Val', -rate, 'name', 'P', 'radius', .1, 'Comp_i', [0 1]);

x  = initState(G,W,100*barsa, [0 1]);
hT = computeTrans(G, rock);
for i=1:nstep
    x  = incompTPFA(x, G, hT, fluid, 'wells', W);
    x  = explicitTransport(x, G, 1/nstep, rock, fluid, 'wells', W);
end

% Set plotting function for later
N = 10;
cval = linspace(0,1,N+1); cval=.5*cval(1:end-1)+.5*cval(2:end);
plotData = @(G,x) contour(reshape(G.cells.centroids(:,1),G.cartDims), ...
                          reshape(G.cells.centroids(:,2),G.cartDims), ...
                          reshape(x,G.cartDims), cval, '-k', 'LineWidth',1.5);
colormap(flipud(.5*jet(N)+.5*ones(N,3)));

% -------------------------------------------------------------------------
% Rotated quater-five spot
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
    xr  = explicitTransport(xr, Gr, 1/nstep, rockr, fluid, 'wells', Wr);
end

% -------------------------------------------------------------------------
% Plot and compare the two solutions
remapAndPlot(Gr,xr.s(:,1),domain([1 2]),cval);
hold on; plotData(G,x.s(:,1)); hold off;
axis equal tight; axis([0 domain(1) 0 domain(2)]);
end