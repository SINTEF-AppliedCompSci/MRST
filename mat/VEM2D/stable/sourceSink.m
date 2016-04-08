clc; clear all; close all;

% addpath('../');
% addpath('../../');
% addpath('~/Documents/master/pebiGridding/voronoi2D/')
% run('../../startup.m');


run('~/NTNU/5/master/project-mechanics-fractures/mystartup.m')
addpath('~/NTNU/5/master/flowSimulationUsingMRST/pebi/')

n = 5;

% x     = linspace(0.2, 0.8, 10);
% y     = 0.8 - 0.5*x - 0.05* sin(6*pi*x);
% fault = {[x' , y']};
% G = compositePebiGrid(1/20, [1, 1], ...
%                        'faultLines', fault, 'faultGridFactor', 1/sqrt(2));

x = linspace(.1,.8,10);
y = 1-x;
fault = {[x' , y']};
G = pebiGrid(1/10, [1, 1], ...
                       'faultLines', fault);
                   
                   
faultFaces = 1:G.faces.num;
faultFaces = faultFaces(G.faces.tag);
G = mrstGridWithFullMappings(G);
G = sortEdges(G);
G = VEM2D_makeInternalBoundary(G, faultFaces);
f = zeros(G.cells.num,1);
G = computeVEM2DGeometry(G,f,1,1);
find(min(abs(bsxfun(@minus,G.cells.centroids,[.2,.2]))));

sourceCoords = [.2,.2];
source = sum(bsxfun(@minus, G.cells.centroids, sourceCoords).^2,2);
source = find(source == min(source));

sinkCoords = [.8,.8];
sink = sum(bsxfun(@minus, G.cells.centroids, sinkCoords).^2,2);
sink = find(sink == min(sink));

f = @(X) zeros(size(X,1),1);

Q = 10;
src = addSource([], source, Q);
src = addSource(src, sink, -Q);

boundaryEdges = find(G.faces.neighbors(:,1) == 0 | G.faces.neighbors(:,2) == 0);
isExternal = G.faces.centroids(boundaryEdges,1) == 0 | ...
        G.faces.centroids(boundaryEdges,1) == 1 | ...
        G.faces.centroids(boundaryEdges,2) == 0 | ...
        G.faces.centroids(boundaryEdges,2) == 1;
isInternal = ~isExternal;

bc = VEM_addBC(G, [], boundaryEdges(isExternal), 'pressure', 0);
bc = VEM_addBC(G, bc, boundaryEdges(isInternal), 'flux', 0);
            
sol1 = VEM2D_v3(G,f,1,bc,'src', src);
sol1 = cellAverages(G,sol1);

% plotVEM(G, sol1.nodeValues, 'dof')

sol2 = VEM2D_v3(G,f,2,bc, 'src', src);

% plotVEM(G,[sol2.nodeValues; sol2.edgeValues; sol2.cellMoments], '')

figure;
plotCellData(G, sol1.cellAverages);
colorbar;
figure;
plotCellData(G, sol2.cellMoments)
colorbar;

% %% Set simulation parameters
% T      = 120*second();    % End time
% dT     = T/120;           % Time steps
% dTplot = 1:1:T;           % Plot at these time steps
% 
% %% Generate grid
% 
% %% Set fluid and rock properties
% gravity reset off 
% 
% fluid = initSimpleFluid('mu' , [   1,   5]*centi*poise     , ...
%                         'rho', [1000, 700]*kilogram/meter^3, ...
%                         'n'  , [   2,   2]);
% 
% rock.poro = ones(G.cells.num,1)*0.15;
% rock.perm = ones([G.cells.num,1])*100*milli*darcy;
% 
% %% add Sources
% srcCells = find(G.cells.tag);
% pv = sum(poreVolume(G, rock));
% 
%             
% %% Solve
% state = initState(G, [], 0, [0.0,1]);
% trans = computeTrans(G,rock);
% % state = incompTPFA(state, G, trans, fluid, 'src', src);
% sol1 = VEM2D_v3(G,f,1,bc,'src', src);
% sol1 = cellAverages(G,sol1);
% state.pressure = sol1.cellAverages;
% % Prepare plotting of saturations
% clf;
% hold on
% plotGrid(G, 'FaceColor', 'none', 'EdgeAlpha', 0.1);
% axis off equal, view([-120,30]), colormap(flipud(jet))
% 
% colorbar; hs = []; ha=[]; zoom(1.3);
% 
% % Start the main loop
% t  = 0;  plotNo = 1;
% while t < T,
%    state = implicitTransport(state, G, dT, rock, fluid, 'src', src);
% 
%    % Check for inconsistent saturations
%    assert(max(state.s(:,1)) < 1+eps && min(state.s(:,1)) > -eps);
%    
%    % Update solution of pressure equation.
%    %state = incompMimetic(state, G3D, S, fluid, 'wells', W);
% %    state = incompTPFA(state, G, trans, fluid, 'src', src);
%     sol1 = VEM2D_v3(G,f,1,bc,'src', src);
%     sol1 = cellAverages(G,sol1);
%     state.pressure = sol1.cellAverages;
%    
%    % Increase time and continue if we do not want to plot saturations
%    t = t + dT;   
%    if ( t + dT <= dTplot(plotNo)), continue, end
%    delete([hs, ha])
%    hs = plotCellData(G, state.s(:,1), find(state.s(:,1) >= 0.0));
%    ha = annotation('textbox', [0.1714 0.8214 0.5000 0.1000], 'LineStyle', 'none', ...
%                    'String', ['Water saturation at ', ...
%                               num2str(convertTo(t,second)), ' s']);
%    fig = gcf();
%    set(findall(fig,'-property','FontSize'),'FontSize',14) 
%    view(0, 90), drawnow, caxis([0 1])
%   
%    plotNo = plotNo+1;
% end