clc; clear; close all;

%% LOAD MODULES

mrstModule add mimetic
mrstModule add mpfa
addpath('../VEM2D/')
addpath('../../../pebiGridding/voronoi2D/')
        
%%  LOAD GRID

grid = 1;

switch grid
    
    case 1
    
        load('singularityGrid.mat')
        xMax = 1; yMax = 1;
        C = [xMax,yMax]/2;
        
    case 2
        
        load('singularityGrid2.mat')
        xMax = 5; yMax = 1;
        C = [xMax,yMax]/2;
        
end

G = computeVEM2DGeometry(G);
        
%%  SET BCs
add = 10;
gD = @(X) -log(sqrt(sum(bsxfun(@minus,X,C).^2,2)))/(2*pi) + add;
gN = @(X) -(X(:,2)-C(2))./(2*pi*sum(bsxfun(@minus,X,C).^2,2)) ;

tol = 1e-6;

boundary = any(G.faces.neighbors == 0,2);
bDir     = find(abs(G.faces.centroids(boundary,1))      < tol | ...
                abs(G.faces.centroids(boundary,1)-xMax) < tol);
bNeuS    = find(abs(G.faces.centroids(boundary,2))      < tol);
bNeuN    = find(abs(G.faces.centroids(boundary,2)-yMax) < tol);
                
bc_VEM = VEM2D_addBC([]    , G, bDir , 'pressure', gD         );
bc_VEM = VEM2D_addBC(bc_VEM, G, bNeuS, 'flux'    , gN         );
bc_VEM = VEM2D_addBC(bc_VEM, G, bNeuN, 'flux'    , @(X) -gN(X));

bc_MRST = addBC([]     , bDir , 'pressure', ...
                                          gD(G.faces.centroids(bDir ,:)) );
bc_MRST = addBC(bc_MRST, bNeuS, 'flux'    , ...
                   -G.faces.areas(bNeuS).*gN(G.faces.centroids(bNeuS,:)) );
bc_MRST = addBC(bc_MRST, bNeuN, 'flux'    , ...
                    G.faces.areas(bNeuN).*gN(G.faces.centroids(bNeuN,:)) );

%% FLUID AND ROCK PROPERTIES

gravity reset off 
fluid = initSingleFluid('mu' , 1, 'rho', 1);
rock.poro = ones(G.cells.num,1);
rock.perm = ones([G.cells.num,1]);

%% INITIALIZE STATE

sInit     = initState(G, [], 0, [0.0,1]);
S         = computeMimeticIP(G, rock, 'Verbose', true);
transTPFA = computeTrans(G,rock);
transMPFA = computeMultiPointTrans(G, rock);

%% SOLVE

sVEM1 = VEM2D(G, 0, bc_VEM, 1, 'cellAverages', true);
sVEM2 = VEM2D(G, 0, bc_VEM, 2                      );

sTPFA = incompTPFA(     sInit, G, transTPFA, fluid, 'bc', bc_MRST);
sMFD  = solveIncompFlow(sInit, G, S        , fluid, 'bc', bc_MRST);
sMPFA  = incompMPFA(    sInit, G, transMPFA, fluid, 'bc', bc_MRST);

sol = [sTPFA.pressure   , ...
       sMPFA.pressure   , ...
       sMFD.pressure    , ...
       sVEM1.cellMoments, ...
       sVEM2.cellMoments] - 10;

%% PLOT

rL = 0:0.001:1;
r = sqrt(sum(bsxfun(@minus,G.cells.centroids,C).^2,2));

plot(rL,gD(XL)-add)
set(gcf, 'DefaultTextInterpreter', 'Latex')
hold on
plot(r, sol(:,1), r, '.', sol(:,1), '.', r,sol(:,3), '.', ...
     r, sol(:,4), '.', r, sol(:,5), '.');
yLim = max(sol)-min(sol);
xLim = max(r)-min(r);
axis([min(r)-xLim*.1, max(r)+xLim*.1, min(sol)-yLim*.1, max(sol)+yLim*.1]);

h = legend('Exact solution', 'TPFA', 'MPFA', 'MFD', '1st order VEM', '2nd order VEM');
set(h, 'interpreter', 'latex');
xlabel('$r$'); ylabel('$p$')
