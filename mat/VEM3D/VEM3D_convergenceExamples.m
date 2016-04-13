%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------

clc; clear; close all;

ex = 1;

switch ex
    case 1
        xMax = 1; yMax = 1; zMax = 1;
        gridType = 'cart';
        neuEx = false;
        f = 0;
        gD = @(X) X(:,1).*X(:,2);
end


nVec = [5,10,20];
nIt = numel(nVec);
errVec = zeros(nIt, 3);

for i = 1:nIt
    
    n = nVec(i);
    if strcmp(gridType,'cart')
        G = cartGrid([n,n,n],[xMax, yMax, zMax]);
    else
        G = unitSquare([n,n,n], [xMax, yMax, zMax]);
    end
    G = computeVEM3DGeometry(G);

    boundaryEdges = find(any(G.faces.neighbors == 0,2));
    if neuEx
        tol = 1e-10;
        isNeu = abs(G.faces.centroids(boundaryEdges,1)) < tol;
    else
        isNeu = false(numel(boundaryEdges),1);
        gN = gD;
    end
    bc = VEM3D_addBC([], boundaryEdges(~isNeu), 'pressure', gD);
    bc = VEM3D_addBC(bc, boundaryEdges(isNeu), 'flux', gN);
    
    [sVEM1, G1] = VEM3D(G,f,bc,1, 'projectors', true);
    [sVEM2, G2] = VEM3D(G,f,bc,2,'projectors', true);
    
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