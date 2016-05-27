clc; clear; close all;

mrstModule add mimetic
mrstModule add mpfa
addpath('../VEM2D/')
addpath('../../../pebiGridding/voronoi2D/')

%% Chose grid
% Chose the type of grid.
% gT = 1      Coarse cartesian
% gT = 2      composite pebi
% gT = 3      fully unstructured grid

xMax = 1; yMax = 1;

C = [xMax/2, yMax/2];
m = 5;
wellLine = {C};                % Set source center

gT = 6;

switch gT
    case 1

        nx = 61;
        G = cartGrid([nx,nx],[xMax,yMax]);
        G = computeGeometry(G);
        w1 = wellLine{1};
        D = pdist2(G.cells.centroids, w1);
        [~, I] = min(D, [], 1);
        G.cells.tag = false(G.cells.num,1);
        G.cells.tag(I') = true(size(I'));
        G = sortEdges(G);
        G = computeVEM2DGeometry(G);

    case 2 

        gridSize = xMax/20;                   % Size of gridcells
        mlqtMax = 2;                            % Set number of reminement levels
        wellGridSize = 0.75/2^mlqtMax;          % Set relative well grid size
        mlqtSizes = 2.0*linspace(gridSize,gridSize*wellGridSize,mlqtMax+1)';
                                                % Size around wells to be refined 
        G = compositePebiGrid(gridSize, [xMax, yMax], 'wellLines', wellLine, ...
                             'wellGridFactor', wellGridSize, ...
                             'mlqtMaxLevel', 2, 'mlqtLevelSteps', mlqtSizes);
        G = sortEdges(G);
        G = computeVEM2DGeometry(G);

    case 3
        
        gridSize = xMax/10;
        wellGridSize = 0.7/2^2;
        epsilon = gridSize*.7;
        G = pebiGrid(gridSize, [xMax, yMax], 'wellLines', wellLine,   ...
                    'wellGridFactor', wellGridSize, 'wellRefinement',true);
        G = sortEdges(G);
        m = 1;
        G = computeVEM2DGeometry(G);

        save('singularityGrid.mat','G');

    case 4

        gridSize = xMax/10;
        wellGridSize = 0.7/2^2;
        epsilon = gridSize*.7;
        G = pebiGrid(gridSize, [xMax, yMax], 'wellLines', wellLine,   ...
                    'wellGridFactor', wellGridSize, 'wellRefinement',true);
        G = sortEdges(G);
        G.nodes.coords(:,1) = G.nodes.coords(:,1)*m;
        xMax = xMax*m;
        C = [xMax/2, yMax/2];
        G = computeVEM2DGeometry(G);
       
        save('singularityGrid2.mat','G');   

    case 5
        
        m = 1;
        load('singularityGrid.mat')
        
    case 6
        
        xMax = xMax*m;
        C = [xMax/2, yMax/2];
        load('singularityGrid2.mat');
       
    
    case 7
        xMax = xMax*1000;
        yMax = xMax;
        C = [xMax/2, yMax/2];
        G = cartGrid([10,10]*4+1,[xMax,yMax]);
        Gvem = computeVEM2DGeometry(G);
        G = computeGeometry(G);
        G.cells.tag(G.cells.centroids(:,1) == C(1) & G.cells.centroids(:,2) ==C(2)) = 1;
        
    otherwise

    error('unknown grid case')

end                  


%%  Set BC
gD = @(X) -1/(2*pi)*log(sqrt(sum(bsxfun(@minus, X, C).^2,2)))+10;
gN = @(X) -1/(2*pi)*(X(:,2)-C(2))./(sum(bsxfun(@minus, X, C).^2,2));
% ggD = @(X) 1/(4*pi*sqrt(bsxfun(@minus,C
% gNn = @(X) -1/(2*pi)*bsxfun(@rdivide,bsxfun(@minus,X,C),(sum(bsxfun(@minus, X, C).^2,2)));

tol = 1e-6;
% % boundaryEdgesDir = find(G.faces.neighbors(:,1) == 0);
% boundary= find(any(G.faces.neighbors == 0,2));
% sign=2*(G.faces.neighbors(boundary,2)==0)-1;
% 
% flux= bsxfun(@times,gNn(G.faces.centroids(boundary,:)).*G.faces.normals(boundary,:),sign)
% sum(sum(flux,2))
%%
boundaryEdgesDir = find(abs(G.faces.centroids(:,1))<tol  | abs(G.faces.centroids(:,1) -xMax)<tol);
boundaryEdgesNeuN = find(abs(G.faces.centroids(:,2)- yMax)<tol);
boundaryEdgesNeuS = find(abs(G.faces.centroids(:,2)) <tol);
bc_VEM = VEM2D_addBC([], G, boundaryEdgesDir, 'pressure', gD);
bc_VEM = VEM2D_addBC(bc_VEM,G,boundaryEdgesNeuN, 'flux', gN);
bc_VEM = VEM2D_addBC(bc_VEM,G,boundaryEdgesNeuS, 'flux', @(X)-gN(X));
bc_MRST = addBC([], boundaryEdgesDir, 'pressure', gD(G.faces.centroids(boundaryEdgesDir,:)));
bc_MRST = addBC(bc_MRST, boundaryEdgesNeuN, 'flux', -2*abs(G.faces.areas(boundaryEdgesNeuN).*gN(G.faces.centroids(boundaryEdgesNeuN,:))));
bc_MRST = addBC(bc_MRST, boundaryEdgesNeuS, 'flux', -2*abs(G.faces.areas(boundaryEdgesNeuS).*gN(G.faces.centroids(boundaryEdgesNeuS,:))));

%% Set fluid and rock properties
gravity reset off 
fluid = initSingleFluid('mu' , 1, 'rho', 1);
rock.poro = ones(G.cells.num,1);
rock.perm = ones([G.cells.num,1]);

%% add Sources
Q = 1;
srcCells = find(G.cells.tag);
src = addSource([],srcCells(1),Q);
src2 = addSource([],srcCells(1),Q);
% W=addWell([],G,rock,srcCells(1),'Type','rate','Val',2*Q,'Radius',0.1/G.cartDims(1));

%% Initialize state
sInit = initState(G, [], 0, [0.0,1]);
S     = computeMimeticIP(G, rock, 'Verbose', true);
trans = computeTrans(G,rock);
transMPFA = computeMultiPointTrans(G, rock);

%% Solve Laplace
sTPFA = incompTPFA(sInit, G, trans, fluid, 'src', src2, 'bc', bc_MRST);
% sTPFA = incompTPFA(sInit, G, trans, fluid, 'wells', W, 'bc', bc_MRST);
sMIM  = solveIncompFlow(sInit, G, S, fluid,'src', src2, 'bc', bc_MRST);
sMPFA  = incompMPFA(sInit, G, transMPFA, fluid, 'src', src2, 'bc',bc_MRST);
sVEM1 = VEM2D(G,0,bc_VEM, 1,'src', src, 'cellAverages', true);
sVEM2 = VEM2D(G,0,bc_VEM, 2,'src', src);


%%  Plot solutions along radial line

dest = '../../tex/thesis/fig/';

r = .01;
n = 20;

xL = linspace(xMax/2*1.05,xMax,n);
yL = yMax/xMax*xL;
XL = [xL', yL'];
rCells = G.cells.centroids(:,1)>xMax/2;
lineCells = zeros(n,1);

for i = 1:n
    dist = sum(bsxfun(@minus,G.cells.centroids,XL(i,:)).^2,2);
    c = find(dist == min(dist));
    lineCells(i) = c(1);
end

XC = G.cells.centroids(lineCells,:);
rC = sqrt(sum(bsxfun(@minus,XC,C).^2,2));
VEM1Vals = sVEM1.cellMoments(lineCells) -10;
VEM2Vals = sVEM2.cellMoments(lineCells) -10;

TPFAVals = sTPFA.pressure(lineCells) -10;
MPFAVals = sMPFA.pressure(lineCells) -10;
MIMVals = sMIM.pressure(lineCells)   -10;

xL = xMax/2:.001:2*xMax;
yL = yMax/xMax*xL;
XL = [xL', yL'];
rL = sqrt(sum(bsxfun(@minus, XL, C).^2,2));

% figure;
% plot(rL,gD(XL)-10)
% set(gcf, 'DefaultTextInterpreter', 'Latex')
% hold on
% plot(rC, TPFAVals, '+', rC, MPFAVals, 'd', rC, MIMVals, 'o', rC, VEM1Vals, '.', rC, VEM2Vals, 'sq');
% h = legend('Exact solution', 'TPFA', 'MPFA', 'MFD', '1st order VEM', '2nd order VEM');
% set(h, 'interpreter', 'latex');
% xlabel('$r$'); ylabel('$p$')
% % set(gca, 'XTick', 0:.1:.7)
% yValMin = min([TPFAVals; MPFAVals; MIMVals; VEM1Vals; VEM2Vals]);
% yValMax = max([TPFAVals; MPFAVals; MIMVals; VEM1Vals; VEM2Vals]);
% yLim = yValMax-yValMin;
% axis([-.1, max(rL), yValMin-yLim*.1, yValMax+yLim*.1]);

%%

r = sqrt(sum(bsxfun(@minus,G.cells.centroids,C).^2,2));
plot(rL,gD(XL)-10)
set(gcf, 'DefaultTextInterpreter', 'Latex')
hold on
plot(r, sTPFA.pressure-10,'.',r, sMPFA.pressure-10, '.',r, sMIM.pressure-10, '.',r, sVEM1.cellMoments-10, '.',r, sVEM2.cellMoments-10, '.')
yValMin = min([sTPFA.pressure; sMPFA.pressure; sMIM.pressure; sVEM1.cellMoments; sVEM2.cellMoments]-10);
yValMax = max([sTPFA.pressure; sMPFA.pressure; sMIM.pressure; sVEM1.cellMoments; sVEM2.cellMoments]-10);
yLim = yValMax-yValMin;
xLim = max(r)-min(r);
axis([min(r)-xLim*.1, max(r)+xLim*.1, yValMin-yLim*.1, yValMax+yLim*.1]);

h = legend('Exact solution', 'TPFA', 'MPFA', 'MFD', '1st order VEM', '2nd order VEM');
set(h, 'interpreter', 'latex');
xlabel('$r$'); ylabel('$p$')

savePdf(gcf, strcat(dest, 'PointSourceComp', num2str(m), '.pdf'));

%%

figure;
plotGrid(G, 'faceAlpha', .2);
set(gcf, 'DefaultTextInterpreter', 'Latex')
% hold on
% % XL = repmat((xMax/2:.001:xMax)',1,2);
% plot(XL(:,1),XL(:,2),'linewidth', 1);
% h1 = plot(XC(:,1), XC(:,2), '.r');
xlabel('$x$'); ylabel('$y$')
set(gca, 'XTick', 0:xMax/2:xMax)
set(gca, 'YTick', 0:yMax/2:yMax)
    
if gT == 3 || gT == 5
    axis equal
    axis([0, 1, 0, 1]);
    savePdf(gcf, strcat(dest, 'PointSourceGrid.pdf'));
end