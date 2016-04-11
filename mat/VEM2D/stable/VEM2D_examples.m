%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------

clc; clear; close all;

ex = 4;
gridType = 'cart';
n = 20;
k = 1;

neuEx = false;
knownSol = true;
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
        f = @(X) pi^2*X(:,1).*sin(pi*X(:,2));
        gD = @(X) X(:,1).*sin(pi*X(:,2));
        gN = @(X) -sin(pi*X(:,2));
end

if strcmp(gridType,'cart')
    G = cartGrid([n,n],[1,1]);
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
    gN = 0;
end
bc = VEM_addBC([], G, boundaryEdges(~isNeu), 'pressure', gD);
bc = VEM_addBC(bc, G, boundaryEdges(isNeu), 'flux', gN);

[sol, G] = VEM2D_v3(G,f,k,bc, 'cellAverages', true);

plotCellData(G, sol.cellMoments)
colorbar;

if knownSol
    l2Err = l2Error(G,sol,gD,k);
    if k == 1
        u = gD(G.nodes.coords);
        err = u - sol.nodeValues;
    elseif k == 2
        u = [gD(G.nodes.coords); gD(G.faces.centroids); polygonInt_v2(G, gD, 7)./G.cells.volumes];
        err = [sol.nodeValues; sol.edgeValues; sol.cellMoments] - u;
    end
    fprintf('l2 error: %d \n\n', sqrt(sum(l2Err)));
    fprintf('2norm error: %d \n\n', norm(err));
end