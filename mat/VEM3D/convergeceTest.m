%--------------------------------------------------------------------------
%   Convergence test.
%--------------------------------------------------------------ØSK-2016021-

clc; clear; close all;
addpath('/home/strene/Documents/master/coop/pebiGridding/voronoi3D')

ex = 3;
switch ex
    case 1
        f  = @(X) pi^2*X(:,1).*sin(pi*X(:,2).*X(:,3)).*(X(:,2).^2 + X(:,3).^2);
        gD = @(X) X(:,1).*sin(pi*X(:,2).*X(:,3));
        grid = 'pebi';
        neu = false;
    case 2
        f  = @(X) pi^2*X(:,1).*sin(pi*X(:,2).*X(:,3)).*(X(:,2).^2 + X(:,3).^2);
        gD = @(X) X(:,1).*sin(pi*X(:,2).*X(:,3));
        grid = 'cart';
        neu = false;
    case 3
        f = @(X) sin(X(:,1));
        gD = @(X) f(X);
        grid = 'pebi';
        neu = false;
end

nVec = [250, 500, 1000, 2000];
nn = numel(nVec);
errVec = zeros(nn, 3);
err = zeros(nn,2);

for i = 1:nn
    
n = nVec(i);

if strcmp(grid, 'cart')
    n = round(n^(1/3));
    G = cartGrid([n,n,n], [1,1,1]);
else
    G = voronoiCube(n,[1,1,1]);
end

%   Compute VEM geometry.
G = computeVEM3DGeometry(G);

%   Define boundary faces.                          
boundaryFaces = find(any(G.faces.neighbors == 0,2));
isNeu = false(numel(boundaryFaces),1);
if neu
    isNeu = boundaryFaces(G.faces.centroids(boundaryFaces,1) == 0);
else
    gN = 0;
end

%   Set boundary conditions.
bc = VEM3D_addBC([], boundaryFaces(~isNeu), 'pressure', gD);
% bc = VEM3D_addBC(bc, boundaryFaces(isNeu), 'flux', gN);

%   Solve.
[sol1, G1] = VEM3D(G, f, bc, 1, 'cellProjectors', true);
[sol2, G2] = VEM3D(G, f, bc, 2, 'cellProjectors', true);

%   L2 error estimate.
l2Err1 = l2Error3D(G1, sol1, gD, 1);
l2Err2 = l2Error3D(G2, sol2, gD, 2);
h = mean(G.cells.diameters);
errVec(i,:) = [h, sqrt(sum(l2Err1)), sqrt(sum(l2Err2))];

end

loglog(errVec(:,1), errVec(:,2), 'sq-', errVec(:,1), errVec(:,3), 'sq-');
s1 = polyfit(log(errVec(:,1)), log(errVec(:,2)), 1);
s2 = polyfit(log(errVec(:,1)), log(errVec(:,3)), 1);
fprintf('Empiric Convergence rate %f for 1st order VEM\n', s1(1));
fprintf('Empiric Convergence rate %f for 2nd order VEM\n\n', s2(1));
% 
% loglog(err(:,1), err(:,2));
% 
% s1 = polyfit(log(err(:,1)), log(err(:,2)), 1);

