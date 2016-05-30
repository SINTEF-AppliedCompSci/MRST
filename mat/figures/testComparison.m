clc; clear; close all;

%% LOAD MODULES

mrstModule add mimetic
mrstModule add mpfa
addpath('../VEM2D/')
addpath('../../../pebiGridding/voronoi2D/')
        
%%  LOAD GRID

grid = 2;

switch grid
    
    case 1
    
        load('singularityGrid.mat')
        xMax = 1; yMax = 1;
        C = [xMax,yMax]/2;
        
    case 2
        
        load('singularityGrid2.mat')
        xMax = 10; yMax = 1;
        C = [xMax,yMax]/2;
        
end

G = computeVEM2DGeometry(G);
        
%%  SET BCs

add = 10;
gD = @(X) -log(sqrt(sum(bsxfun(@minus,X,C).^2,2)))/(2*pi) + add;
gN = @(X) -(X(:,2)-C(2))./(2*pi*sum(bsxfun(@minus,X,C).^2,2)) ;

tol = 1e-6;

bDir     = find(abs(G.faces.centroids(:,1))      < tol | ...
                abs(G.faces.centroids(:,1)-xMax) < tol);
bNeuS    = find(abs(G.faces.centroids(:,2))      < tol);
bNeuN    = find(abs(G.faces.centroids(:,2)-yMax) < tol);
                
bc_VEM = VEM2D_addBC([]    , G, bDir , 'pressure', gD         );
bc_VEM = VEM2D_addBC(bc_VEM, G, bNeuS, 'flux'    , @(X) -gN(X));
bc_VEM = VEM2D_addBC(bc_VEM, G, bNeuN, 'flux'    , gN         );

bc_MRST = addBC([]     , bDir , 'pressure', ...
                                          gD(G.faces.centroids(bDir ,:)) );
bc_MRST = addBC(bc_MRST, bNeuS, 'flux'    , ...
                   -G.faces.areas(bNeuS).*gN(G.faces.centroids(bNeuS,:)) );
bc_MRST = addBC(bc_MRST, bNeuN, 'flux'    , ...
                    G.faces.areas(bNeuN).*gN(G.faces.centroids(bNeuN,:)) );

%% FLUID AND ROCK PROPERTIES

gravity off 
fluid = initSingleFluid('mu' , 1, 'rho', 1);
rock = struct('perm', ones(G.cells.num,1), 'poro', ones(G.cells.num,1));
% rock.poro = ones(G.cells.num,1);
% rock.perm = ones([G.cells.num,1]);

%% ADD SOURCE

Q = 1;
srcCells = find(G.cells.tag);
src = addSource([],srcCells(1),Q);
src2 = addSource([],srcCells(1),1*Q);

%% INITIALIZE STATE

sInit     = initState(G, [], 0);
S         = computeMimeticIP(G, rock, 'Verbose', true);
transTPFA = computeTrans(G,rock);
transMPFA = computeMultiPointTrans(G, rock);

%% SOLVE

sTPFA = incompTPFA(     sInit, G, transTPFA, fluid, 'bc', bc_MRST, 'src', src2);
sMPFA  = incompMPFA(    sInit, G, transMPFA, fluid, 'bc', bc_MRST, 'src', src2);
sMFD  = solveIncompFlow(sInit, G, S        , fluid, 'bc', bc_MRST, 'src', src2);

sVEM1 = VEM2D(G, 0, bc_VEM, 1, 'src', src, 'cellAverages', true);
sVEM2 = VEM2D(G, 0, bc_VEM, 2, 'src', src                      ); 

sol = [sTPFA.pressure   , ...
       sMPFA.pressure   , ...
       sMFD.pressure    , ...
       sVEM1.cellMoments, ...
       sVEM2.cellMoments] - add;

%% PLOT

dest = '../../tex/thesis/fig/';

r = sqrt(sum(bsxfun(@minus,G.cells.centroids,C).^2,2));
xL = linspace(xMax/2,2*xMax,1000);
yL = yMax/xMax*xL;
XL = [xL', yL'];
rL = sqrt(sum(bsxfun(@minus, XL, C).^2,2));

plot(rL,gD(XL)-add)
set(gcf, 'DefaultTextInterpreter', 'Latex')
hold on
plot(r, sol(:,1), '.', r, sol(:,2), '.', r,sol(:,3), '.', ...
     r, sol(:,4), '.', r, sol(:,5), '.');
yMaxL = mean(max(sol)); yMinL = mean(min(sol));
yLim = yMaxL-yMinL;
xLim = max(r)-min(r);
axis([min(r)-xLim*.1, max(r)+xLim*.1, yMinL-yLim*.1, .67]);

h = legend('Exact solution', 'TPFA', 'MPFA', 'MFD', '1st order VEM', '2nd order VEM');
set(h, 'interpreter', 'latex');
xlabel('$r$'); ylabel('$p$')

savePdf(gcf, strcat(dest, 'PointSourceComp', num2str(grid), '.pdf'));

if grid == 1
    figure;
    plotGrid(G, 'facealpha', .2);
    hold on
    plotGrid(G, srcCells, 'facecolor', [0.8500 0.3250 0.0980]);
    set(gcf, 'DefaultTextInterpreter', 'Latex')
    xlabel('$x$'); ylabel('$y$')
    set(gca, 'XTick', 0:xMax/2:xMax)
    set(gca, 'YTick', 0:yMax/2:yMax)
    axis equal
    axis([0, 1, 0, 1]);
    savePdf(gcf, strcat(dest, 'PointSourceGrid.pdf'));
end
