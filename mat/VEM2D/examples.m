%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------

clc; clear; close all;

ex = 3;
gridType = 'pebi';
n = 10;
nx = n; ny = n;
xMax = 1; yMax = 1;
k = 1;

neuEx = false;
knownSol = true;
edgeclr = 'k';
simga = 3;
switch ex
    case 1
        edgeclr = 'none';
        f = 0;
        C = -[.2,.2];
        gD = @(X) -log(1./(sqrt(sum(bsxfun(@minus, X, C).^2,2))));
    case 2
        xMax = 1; yMax = 1;
%         ny = 5*ny;
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
    case 5
        gridType = 'cart';
        nx = 10; ny = 40;
        knownSol = false;
        epsilon = .5*1e-4;
        f = @(X) 1/(2*pi*epsilon)*exp(-((X(:,1)-xMax/2).^2 + (X(:,2)-yMax/2).^2)/(2*epsilon));
        gD = 0;
    case 6
        gridType = 'cart';
        nx = 10; ny = 50;
        k = 1;
        f = @(X) pi^2*X(:,1).*sin(pi*X(:,2));
        gD = @(X) X(:,1).*sin(pi*X(:,2));
        nk   = (k+1)*(k+2)/2;
        NK   = (4 + 4*(k-1) + k*(k-1)/2)*ones(nx*ny,1);
        nker = NK - nk;
        hx = xMax/nx; hy = yMax/ny;
end

if strcmp(gridType,'cart')
    G = cartGrid([nx,ny],[xMax, yMax]);
else
    G = unitSquare([nx,ny], [xMax, yMax]);
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
bc = VEM2D_addBC([], G, boundaryEdges(~isNeu), 'pressure', gD);
bc = VEM2D_addBC(bc, G, boundaryEdges(isNeu), 'flux', gN);

[sol, G] = VEM2D(G,f,k,bc, 'cellAverages', true);

figure;
plotVEM2D(G,sol,k)



if knownSol
    l2Err = l2Error(G,sol,gD,k);
    if k == 1
        u = gD(G.nodes.coords);
        err = u - sol.nodeValues;
    elseif k == 2
        u = [gD(G.nodes.coords); gD(G.faces.centroids); polygonInt(G, 1:G.cells.num, gD, 7)./G.cells.volumes];
        err = [sol.nodeValues; sol.edgeValues; sol.cellMoments] - u;
    end
    fprintf('l2 error: %d \n\n', sqrt(sum(l2Err)));
    fprintf('2norm error: %d \n\n', norm(err));
end