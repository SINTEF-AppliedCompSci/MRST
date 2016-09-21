clc; clear all; close all;

%   1) Validation of consistency for first order.
%   2) Validation of consistency for second order.
%   3) Point source problem, both orders.

i = 1;

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
    
        G = unitSquare([5,10],[1,1]);
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
        
    case 3
        
%         G = unitSquare([10,10],[1,1]);
        G = cartGrid([10,10],[1,1]);
        G = computeVEMGeometry(G);
        
        rock.perm = ones(G.cells.num,1);
        fluid = initSingleFluid('mu', 1, 'rho', 1);
        
        
        Q = 100;
        d = sum(bsxfun(@minus, G.cells.centroids, .5*[1,1]).^2,2);
        srcCell = find(d == min(d));
        srcCell = srcCell(1);
        src = addSource([], srcCell, Q);
        xSrc = G.cells.centroids(srcCell, :);
        
        gD = @(x) -log(sqrt(sum(bsxfun(@minus, x, xSrc).^2,2)))*Q/(2*pi);
        gN = @(x) -(x(:,2)-xSrc(2))./(sum(bsxfun(@minus, x, xSrc).^2,2))*Q/(2*pi);
        tol = 1e-6;
        f = boundaryFaces(G);
        isDir  = abs(G.faces.centroids(f,1)) < tol | abs(G.faces.centroids(f,1)-1) < tol;
        isNeuS = abs(G.faces.centroids(f,2)) < tol;
        isNeuN = abs(G.faces.centroids(f,2)-1) < tol;
        
%         bc = addBCFunc([], f(isDir) , 'pressure', gD);
%         bc = addBCFunc(bc, f(isNeuN), 'flux'    , gN);
%         bc = addBCFunc(bc, f(isNeuS), 'flux'    , @(x) -gN(x));
        
        bc = addBC([], f(isDir) , 'pressure', gD(G.faces.centroids(f(isDir),:)));
        bc = addBC(bc, f(isNeuN), 'flux'    , gN(G.faces.centroids(f(isNeuN),:)).*G.faces.areas(f(isNeuN)) );
        bc = addBC(bc, f(isNeuS), 'flux'    ,-gN(G.faces.centroids(f(isNeuS),:)).*G.faces.areas(f(isNeuS)) );

        
        state = initState(G, [], 0);
        k = 1;
        tic;
        S = computeVirtualIP(G, rock, k);
        state = incompVEM(state, G, S, fluid, 'bc', bc, 'src', src);
        toc
        
%         fprintf('\nError: %.2d\n\n', norm(state.nodePressure -  gD(G.nodes.coords))/norm(gD(G.nodes.coords)));
        
        gr = @(r) -log(r)*Q/(2*pi);
        rr = 0:.001:.8;
        rN = sqrt(sum(bsxfun(@minus, G.nodes.coords, xSrc).^2,2));
        
        if k == 1
            plot(rN, state.nodePressure, '.', rr, gr(rr), '-.');
        else
            rE = sqrt(sum(bsxfun(@minus, G.faces.centroids, xSrc).^2,2));
            rC = sqrt(sum(bsxfun(@minus, G.cells.centroids, xSrc).^2,2));
            plot(rN, state.nodePressure, '.', rE, state.facePressure, '+', rC, state.cellPressure, 'd', rr, gr(rr), '-.');
        end
        
    case 4
        
        
        n = 21;
        G = cartGrid([n,n,n], [1,1,1]);
        
%         G = voronoiCube(2000, [1,1,1000]);
        
%         G = computeVEM3DGeometry(G);
        G = computeVEMGeometry(G);
        
        k = 1;
        
        mu = 1; rho = 1;

        rock.perm = ones(G.cells.num,1);
        fluid = initSingleFluid('mu', mu, 'rho', 1);
        state = initState(G, [], 0);
        
        tol = 1e-6;
        f = boundaryFaces(G);
        
%         gD = @(x) x(:,1).^2+x(:,2).^2-2*x(:,3).^2;
        gD = @(x) x(:,1).^3;
        src = addSource([], 1:G.cells.num, -6*G.cells.volumes.*G.cells.centroids(:,1));

        bc = addBC([], f, 'pressure', gD(G.faces.centroids(f,:)));
            
%         bc = addBCFunc([], f, 'pressure', gD);
%         bc = addBCFunc([], f, 'flux', 0);
        tic;
        S = computeVirtualIP(G, rock, k);
        state = incompVEM(state, G, S, fluid, 'bc', bc, 'matrixOutput', true, 'src', src);
        state = calculateCellPressure(state, G, S);
        toc
        
        fprintf('\nError: %.2f %%\n\n', ...
                norm(state.nodePressure -  gD(G.nodes.coords))/norm(gD(G.nodes.coords))*100);
            
        plotCellData(G, state.cellPressure)
        
    case 5
        
        n = 1;
%         G = cartGrid([n,n,n],[1,1,1]);
%         G = computeVEMGeometry(G);
% %         G.nodes.coords(:,1) = G.nodes.coords(:,1) + G.nodes.coords(:,3);
%         G.nodes.coords(:,1) = G.nodes.coords(:,1).^2;
%         G.nodes.coords(:,2) = G.nodes.coords(:,2).^2;
%         G.nodes.coords(:,3) = G.nodes.coords(:,3).^2;
        xx = linspace(0,1,3); yy = linspace(0,1,3); zz = linspace(0, 1, 3);
        [x,y,z] = meshgrid(xx, yy, zz);
        p = [x(:), y(:), z(:)];
        G = tetrahedralGrid(p);
        G = voronoiCube(100, [1,1,1]);
%         
%         p = [0 0 0; 1 0 0; 0 1 0; 0 0 1];
%         G = tetrahedralGrid(p);
        G = computeVEMGeometry(G);
        
        k = 2;

        rock.perm = ones(G.cells.num,1);
        fluid = initSingleFluid('mu', 1, 'rho', 1);
        state = initState(G, [], 0);

        f = boundaryFaces(G);
        
        gD = @(x) x(:,1).^2 - x(:,3).^2;

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
            
        plotCellData(G, state.cellPressure)
            
    case 6


        k = 1;
        
        n = 21;
        xMax = 1; yMax = 1; zMax = 1;
        G = cartGrid([n,n,n], [xMax, yMax, zMax]);
%         G = voronoiCube(1000, [xMax, yMax, zMax]);        
        G = computeVEMGeometry(G);
        
        Q = 1;
        d = sum(bsxfun(@minus, G.cells.centroids, .5*[xMax, yMax, zMax]).^2,2);
        srcCell = find(d == min(d));
        srcCell = srcCell(1);
        src = addSource([], srcCell, Q);
        xSrc = G.cells.centroids(srcCell, :);
        
        rock.perm = ones(G.cells.num,1);
        mu = 1; rho = 1;
        fluid = initSingleFluid('mu', mu, 'rho', rho);
        state = initState(G, [], 0);
        
        [bf, cc] = boundaryFaces(G);
        bn = G.faces.nodes(mcolon(G.faces.nodePos(bf), G.faces.nodePos(bf+1)-1));
        gD = @(x) Q./(4*pi*sqrt(sum(bsxfun(@minus, x, xSrc).^2,2)));
        bc = addBC([], bf, 'pressure', gD(G.faces.centroids(bf,:)));
        bcVEM = addBCFunc([], bf, 'pressure', gD);
        
        tic;
        SVEM = computeVirtualIP(G, rock, k,  'innerProduct', 'ip_simple');
        toc
        stateVEM = incompVEM(state, G, SVEM, fluid, 'src', src, 'bc', bcVEM);
        stateVEM = calculateCellPressure(stateVEM, G, SVEM);
        
        SMFD = computeMimeticIP(G, rock);
        stateMFD = incompMimetic(state, G, SMFD, fluid, 'src', src, 'bc', bc);
        
        %%
        
        figure;
        rN = sqrt(sum(bsxfun(@minus, G.nodes.coords, xSrc).^2,2));
        rC = sqrt(sum(bsxfun(@minus, G.cells.centroids, xSrc).^2,2));
        gr = @(r) Q./(4*pi*r);
        rr = linspace(0, max(rC), 1000);
        
        plot(rN, stateVEM.nodePressure, '.', rC, stateVEM.cellPressure, 'o', rC, stateMFD.pressure, 'd', rr, gr(rr), '-.');
        ylim([0 Q]);
        
        pc = find(G.cells.centroids(:,3) > .5);
        figure;
        plotCellData(G, stateVEM.cellPressure(pc), pc)
        axis equal
                
        tr = true(G.cells.num,1);
        tr(rC < 1*mean(G.cells.diameters)^2) = false;
        err = norm(stateVEM.cellPressure(tr) - gD(G.cells.centroids(tr,:)))/norm(gD(G.cells.centroids(tr,:)));
        fprintf('Error: %.2f %%\n', err*100);

        
        
%         figure;
%         plot(rN(bn), stateVEM.nodePressure(bn), '.', rN(bn), gD(G.nodes.coords(bn,:)), 'sq')
        
end

