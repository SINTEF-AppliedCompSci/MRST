%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------

clc; clear; close all;

ex = 1;
gridType = 'cart';
plotSol = false;

neuEx = false;
switch ex
    case 1
        f = @(X) zeros(size(X,1),1);
        C = -[.5,.5];
        gD = @(X) -log(1./(sqrt(sum(bsxfun(@minus, X, C).^2,2))));
    case 2
        f = @(X) sin(X(:,1));
        gD = @(X) sin(X(:,1));
    case 3
        neuEx = true;
        f = @(X) 24*exp(X(:,1)).*cos(5*X(:,2));
        gD = @(X) exp(X(:,1)).*cos(5*X(:,2));
        gN = @(X) -gD(X);
    case 4
        neuEx = true;
        f = @(X) -6*X(:,1);
        gD = @(X) X(:,1).^3;
        gN = @(X) -3*X(:,1).^2;
    case 5
        neuEx = true;
        f = @(X) pi^2*X(:,1).*sin(pi*X(:,2));
        gD = @(X) X(:,1).*sin(pi*X(:,2));
        gN = @(X) -sin(pi*X(:,2));
end


nVec = [10, 20, 40, 80];
nIt = numel(nVec);
errVec = zeros(nIt, 3);

for i = 1:nIt
    
    n = nVec(i);
    if strcmp(gridType,'cart')
        G = cartGrid([n,n],[1,1000]);
    else
        G = unitSquare([n,n]);
    end
    G = sortEdges(G);
    G = computeVEM2DGeometry(G);

    boundaryEdges = find(any(G.faces.neighbors == 0,2));
    if neuEx
        tol = 1e-10;
        isNeu = abs(G.faces.centroids(boundaryEdges,1)) < tol;
    else
        isNeu = false(numel(boundaryEdges),1);
        gN = gD;
    end
    bc = VEM_addBC([], G, boundaryEdges(~isNeu), 'pressure', gD);
    bc = VEM_addBC(bc, G, boundaryEdges(isNeu), 'flux', gN);
    
    [sVEM1, G1] = VEM2D_v3(G,f,1,bc, 'projectors', true);
    [sVEM2, G2] = VEM2D_v3(G,f,2,bc, 'projectors', true);
    
    h = mean(G.cells.diameters);
    area = sqrt(sum(G.cells.volumes.^2));
    nK = G.cells.num;
        
    l2Err1 = l2Error(G1, sVEM1, gD, 1);
    l2Err2 = l2Error(G2, sVEM2, gD, 2);
    errVec(i,:) = [h, sqrt(sum(l2Err1)), sqrt(sum(l2Err2))];
    
    if plotSol
        subplot(1,2,1)
        plotCellData(G, sVEM1.cellMoments);
        axis equal;
        subplot(1,2,2);
        plotCellData(G, sVEM2.cellMoments);
        axis equal;
        pause;
    end
end

loglog(errVec(:,1), errVec(:,2), errVec(:,1), errVec(:,3));
legend('VEM 1st order', 'VEM 2nd order');
p1 = polyfit(log(errVec(:,1)), log(errVec(:,2)),1);
p2 = polyfit(log(errVec(:,1)), log(errVec(:,3)),1);

p1(1)
p2(1)