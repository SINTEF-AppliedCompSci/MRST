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


nVec = [3,6,12];
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
    G1 = VEM3D_faceProjectors(G,1);
    G2 = VEM3D_faceProjectors(G,2);

    boundaryEdges = find(any(G1.faces.neighbors == 0,2));
    if neuEx
        tol = 1e-10;
        isNeu = abs(G1.faces.centroids(boundaryEdges,1)) < tol;
    else
        isNeu = false(numel(boundaryEdges),1);
        gN = gD;
    end
    bc = VEM3D_addBC([], boundaryEdges(~isNeu), 'pressure', gD);
    bc = VEM3D_addBC(bc, boundaryEdges(isNeu), 'flux', gN);
    
    [sVEM1, G1] = VEM3D(G1,f,bc,1, 'cellProjectors', true);
    [sVEM2, G2] = VEM3D(G2,f,bc,2, 'cellProjectors', true);
    
    h = mean(G.cells.diameters);
    area = sqrt(sum(G.cells.volumes.^2));
    nK = G.cells.num;
        
    l2Err1 = l2Error3D(G1, sVEM1, gD, 1);
    l2Err2 = l2Error3D(G2, sVEM2, gD, 2);
    errVec(i,:) = [h, sqrt(sum(l2Err1)), sqrt(sum(l2Err2))];

end

