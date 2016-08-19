clc; clear; close all;

mrstModule add mimetic
mrstModule add mpfa
addpath('../VEM2D/')

tol= 1e-6;

i = 3;
switch i
    case 1

        load('unstructPebi.mat');

        k = 2;

        rock.perm = ones(G.cells.num,1);
        fluid = initSingleFluid('mu', 1, 'rho', 1);

        C = [.5, .5];
        gD = @(X) -1/(2*pi)*log(sqrt(sum(bsxfun(@minus, X, C).^2,2)))+10;
        gN = @(X) -1/(2*pi)*(X(:,2)-C(2))./(sum(bsxfun(@minus, X, C).^2,2));

        tol = 1e-6;
        f = boundaryFaces(G);
        isNeuN = abs(G.faces.centroids(f,2) - 1)<tol;
        isNeuS = abs(G.faces.centroids(f,2)) <tol;
        fNeuN = f(isNeuN);
        fNeuS = f(isNeuS);
        fDir = f(~any([isNeuN, isNeuS],2));

        bcFunc = addBCFunc([], fDir, 'pressure', gD);
        bcFunc = addBCFunc(bcFunc, fNeuN, 'flux', gN);
        bcFunc = addBCFunc(bcFunc, fNeuS, 'flux', @(X) -gN(X));

        bc = func2mrstBC(bcFunc, G);

        c = 1:G.cells.num;
        src = addSource([], c(G.cells.tag), 1);
    
    case 2
        
        k = 2;
        src = [];
        
        G = unitSquare([20,20],[1,1]);
        G = computeVEM2DGeometry(G);
        sortEdges(G);
        rock.perm = repmat([10, 0, 1]*1e-12, G.cells.num,1);
        mu = 100; rho = 1;
        fluid = initSingleFluid('mu', mu, 'rho', rho);
        
        h = .7;
        f = boundaryFaces(G);
        isIn = G.faces.centroids(f,2) < h & ...
               abs(G.faces.centroids(f,1)) < tol;
        isOut = G.faces.centroids(f,2) > 1-h & ...
                abs(G.faces.centroids(f,1)-1) < tol;
        bc = addBC([], f(isIn), 'pressure', 100);
        bc = addBC(bc, f(isOut), 'pressure', -100);
        bc = addBC(bc, f(~any([isIn, isOut], 2)), 'flux', 0);
         
%         bcFunc = addBCFunc([], f(isIn), 'flux', 100);
%         bcFunc = addBCFunc(bcFunc, f(isOut), 'flux', -100);
%         bcFunc = addBCFunc(bcFunc, f(~any([isIn, isOut], 2)), 'flux', 0);
        
        bcFunc = bc;
        
        sum(bc.value)
        
        bFunc = bc;
        
        figure;
        plotGrid(G);
        hold on;
        
        plotBC(bc,G)
        
    case 3
        
        k = 2;
        src = [];
        
        G = unitSquare([30,30],[1,1]);
        G = computeVEM2DGeometry(G);
        sortEdges(G);
        rock.perm = [rand(G.cells.num,1)*1000,rand(G.cells.num,1)]*1e-12;
        fluid = initSingleFluid('mu', 1, 'rho', 1);
        
        h = .3;
        f = boundaryFaces(G);
        isIn = G.faces.centroids(f,2) < h & ...
               abs(G.faces.centroids(f,1)) < tol;
        isOut = G.faces.centroids(f,2) > 1-h & ...
                abs(G.faces.centroids(f,1)-1) < tol;
        bc = addBC([], f(isIn), 'flux', 10*G.faces.areas(f(isIn)));
        bc = addBC(bc, f(isOut), 'flux', -10*G.faces.areas(f(isOut)));
        bc = addBC(bc, f(~any([isIn, isOut], 2)), 'pressure', 0);
        
        bcFunc = addBCFunc([], f(isIn), 'flux', 10);
        bcFunc = addBCFunc(bcFunc, f(isOut), 'flux', -10);
        bcFunc = addBCFunc(bcFunc, f(~any([isIn, isOut], 2)), 'pressure', 0);
        
        figure;
        plotGrid(G);
        hold on;
        
        plotBC(bc,G)
        
end

figure;


state = initState(G, [], 0);

tic;
S = computeVirtualIP(G, rock, k);
stateVEM = incompVEM(state, G, S, fluid, 'bc', bcFunc, 'src', src);
toc

subplot(2,2,1)
plotCellData(G, stateVEM.cellPressure);
axis([0,1,0,1]);
colorbar;

S = computeMimeticIP(G, rock);
stateMFD = incompMimetic(state, G, S, fluid, 'bc', bc, 'src', src);

subplot(2,2,2)
plotCellData(G, stateMFD.pressure);
axis([0,1,0,1]);
colorbar;

T = computeTrans(G, rock);
stateTPFA = incompTPFA(state, G, T, fluid, 'bc', bc, 'src', src);

subplot(2,2,3)
plotCellData(G, stateTPFA.pressure);
axis([0,1,0,1]);
colorbar;

T = computeMultiPointTrans(G, rock);
stateMPFA = incompMPFA(state, G, T, fluid, 'bc', bc, 'src', src);

subplot(2,2,4)
plotCellData(G, stateMPFA.pressure);
axis([0,1,0,1]);
colorbar;

fprintf('\nDifference\n-----------------\n')
fprintf('MFD \t %.2f%%\n', 100*norm(stateVEM.cellPressure-stateMFD.pressure)/norm(stateMFD.pressure))
fprintf('MPFA \t %.2f%%\n', 100*norm(stateVEM.cellPressure-stateMPFA.pressure)/norm(stateMPFA.pressure))
fprintf('TPFA \t %.2f%%\n\n', 100*norm(stateVEM.cellPressure-stateTPFA.pressure)/norm(stateMPFA.pressure))