clc; clear; close all;

%   1) Validation of consistency for first order.
%   2) Validation of consistency for second order.
%   3) Point source problem, both orders.

i = 5;

switch i
    case 1

        tol= 1e-6;
        G = unitSquare([10,10],[1,1]);
        G = sortEdges(G);
        G = computeVEMGeometry(G);

        state = initState(G, [], 0);

        k = 1;
        
        K = rand(1,3)*1e-12;
        rock.perm = repmat(K, G.cells.num,1);
        mu = 10; rho = 1;
        fluid = initSingleFluid('mu', mu, 'rho', rho);

        tol = 1e-6;
        f = boundaryFaces(G);
        isNeu = abs(G.faces.centroids(f,1))<tol;
        
        gD = @(X) sum(X,2);
        bc = addBCFunc([], f(isNeu), 'flux', -(K(1) + K(2))/mu);
        bc = addBCFunc(bc, f(~isNeu), 'pressure', gD);
        
        tic;
        S = computeVirtualIP(G, rock, k);
        state = incompVEM(state, G, S, fluid, 'bc', bc);
        toc
        
        fprintf('\nError: %.2d\n\n', norm(state.nodePressure -  gD(G.nodes.coords)));
        
        plotVEM2D(G, state, k);

    case 2
    
        tol= 1e-6;
        G = unitSquare([10,10],[1,1]);
        G = sortEdges(G);
        G = computeVEMGeometry(G);

        state = initState(G, [], 0);
        
        k = 2;
        
        K = rand(1,3)*1e-12;
        
        mu = 100; rho = 1;
        gD = @(X) -(K(2) + K(3))/K(1)*X(:,1).^2 + X(:,1).*X(:,2) + X(:,2).^2;
        gN = @(X) -((2*X(:,2)-X(:,1))*K(2) - 2*K(3)*X(:,1) + K(1)*X(:,2))/mu;
        
        
        rock.perm = repmat(K, G.cells.num,1);
        fluid = initSingleFluid('mu', mu, 'rho', 1);

        tol = 1e-6;
        f = boundaryFaces(G);
        isNeu = abs(G.faces.centroids(f,1))<tol;
        
        bc = addBCFunc([], f(isNeu), 'flux', gN);
        bc = addBCFunc(bc, f(~isNeu), 'pressure', gD);
        
        tic;
        S = computeVirtualIP(G, rock, k);
        state = incompVEM(state, G, S, fluid, 'bc', bc);
        toc
        
        fprintf('\nError: %.2d\n\n', norm(state.nodePressure -  gD(G.nodes.coords))/norm(gD(G.nodes.coords)));
        
        plotVEM2D(G, state, k);
        
    case 3
        
        
        n = 10;
        G = cartGrid([n,n,n], [1,1,10000]);
        
%         G = voronoiCube(2000, [1,1,1000]);
        
%         G = computeVEM3DGeometry(G);
        G = computeVEMGeometry(G);
        
        k = 1;
        
        mu = 100; rho = 1;

        rock.perm = ones(G.cells.num,1);
        fluid = initSingleFluid('mu', mu, 'rho', 1);
        state = initState(G, [], 0);
        
        tol = 1e-6;
        f = boundaryFaces(G);
        
        gD = @(x) x(:,3);

        bc = addBCFunc([], f, 'pressure', gD);
%         bc = addBCFunc([], f, 'flux', 0);
        tic;
        S = computeVirtualIP(G, rock, k);
        state = incompVEM(state, G, S, fluid, 'bc', bc);
        toc
        
        fprintf('\nError: %.2d\n\n', ...
                norm(state.nodePressure -  gD(G.nodes.coords)));
        
    case 4
        
        n = 5;
        G = cartGrid([n,n,n],[1,1,1]);
        
        G = voronoiCube(5, [1,1,1]);
        
%         G = computeVEM3DGeometry(G);
        
        G = computeVEMGeometry(G);
        
        k = 2;
        
        mu = 1; rho = 1;

        rock.perm = ones(G.cells.num,1);
        fluid = initSingleFluid('mu', mu, 'rho', 1);
        state = initState(G, [], 0);

        f = boundaryFaces(G);
        
        gD = @(x) x(:,2).^2 - x(:,3).^2;

        bc = addBCFunc([], f, 'pressure', gD);
        tic;
        S = computeVirtualIP(G, rock, k);
        toc
        state = incompVEM(state, G, S, fluid, 'bc', bc);
        
        p = [gD(G.nodes.coords); gD(G.edges.centroids); ...
             polygonInt3D(G, 1:G.faces.num, gD, 2)./G.faces.areas; ...
             polyhedronInt(G, 1:G.cells.num, gD, 2)./G.cells.volumes];
        P = [state.nodePressure; state.edgePressure; state.facePressure; state.cellPressure];
         
        fprintf('\nError: %.2d\n\n', ...
                norm(p-P));
            
    case 5
        
        n = 20;
        G = cartGrid([n,n,n], [1,1,1]);
%         G = voronoiCube(10, [1,1,1]);
        G = computeVEMGeometry(G);
        
        k = 1;
        
        mu = 1; rho = 1;
        rock.perm = ones(G.cells.num,1);
        fluid = initSingleFluid('mu', mu, 'rho', rho);
        state = initState(G, [], 0);

        d = sum(bsxfun(@minus, G.cells.centroids, .25*[1,1,1]).^2,2);
        srcCell = find(d == min(d));
        srcCell = srcCell(1);
        src = addSource([], srcCell, 10000);
        
        d = sum(bsxfun(@minus, G.cells.centroids, .75*[1,1,1]).^2,2);
        srcCell = find(d == min(d));
        srcCell = srcCell(1);
        src = addSource(src, srcCell, -10000);
        
        b = boundaryFaces(G);
%         bc = addBC([], b, 'pressure', 0);
%         bc = [];
                
        tic;
        S = computeVirtualIP(G, rock, k);
        toc
        stateVEM = incompVEM(state, G, S, fluid, 'src', src);
        stateVEM = calculateCellPressure(stateVEM, G, S);
        
        
%         figure;
%         r = sqrt(sum(bsxfun(@minus, G.cells.centroids, [0.5, 0.5, 0.5]).^2,2));
%         plot(r, stateVEM.cellPressure, '.');
%         figure;
%         r = sqrt(sum(bsxfun(@minus, G.nodes.coords, [0.5, 0.5, 0.5]).^2,2));
%         plot(r, stateVEM.nodePressure, '.');
        
        
        S = computeMimeticIP(G, rock);
        stateMFD = incompMimetic(state, G, S, fluid, 'src', src);
        
        fprintf('\nError: %.2f\n\n', norm(stateVEM.cellPressure - stateMFD.pressure)/norm(stateMFD.pressure));
        
%         
end

