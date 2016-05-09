clc; clear; close all;

addpath('../')
addpath('../../../pebiGridding/voronoi3D/')

ex = 4;
switch ex
    case 1
        %   Specify problem
        f  = @(X) pi^2*X(:,1).*sin(pi*X(:,2).*X(:,3)).*(X(:,2).^2 + X(:,3).^2);
        gD = @(X) X(:,1).*sin(pi*X(:,2).*X(:,3));
        neu = false;
        %   Specify grid
        grid = 'pebi';
        xMax = 1; yMax = 1; zMax = 1;
        n = 200;
        %   Method order
        k = 1;
    case 2
        f = @(X) -4*ones(size(X,1),1);
        gD = @(X) X(:,1).^2 + X(:,2).*X(:,3)*10 + X(:,3).^2;
        neu = false;
        grid = 'pebi';
        xMax = 1; yMax = 1; zMax = 1;
        nx = 5; ny = 5; nz = 5;
        n = nx*ny*nz;
        %   Method order
        k = 2;
    case 3
        f = @(X) 0*ones(size(X,1),1);
        gD = @(X) X(:,1) + 100*X(:,3) + X(:,2)/6;
        neu = false;
        grid = 'pebi';
        xMax = 1; yMax = 1; zMax = 1;
        nx = 5; ny = 5; nz = 5;
        n = nx*ny*nz;
        %   Method order
        k = 1;
    case 4
        f = @(X) -4*ones(size(X,1),1);
        gD = @(X) X(:,1).^2 + X(:,2).*X(:,3)*10 + X(:,3).^2;
        gN = @(X) 10*X(:,2) + 2*X(:,3);
        neu = true;
        grid = 'pebi';
        xMax = 1; yMax = 1; zMax = 1;
        nx = 5; ny = 5; nz = 5;
        n = nx*ny*nz;
        %   Method order
        k = 2;
end

if strcmp(grid, 'cart')
    G = cartGrid([nx,ny,nz], [xMax, yMax, zMax]);
else
    G = voronoiCube(n, [xMax, yMax, zMax]);
end

G = computeVEM3DGeometry(G);
                           
boundaryFaces = find(any(G.faces.neighbors == 0, 2));
isNeu = false(numel(boundaryFaces),1);
if neu
    isNeu = G.faces.centroids(boundaryFaces,3) == 1;
else
    gN = 0;
end

%   Set boundary conditions.
bc = VEM3D_addBC([], boundaryFaces(~isNeu), 'pressure', gD);
bc = VEM3D_addBC(bc, boundaryFaces(isNeu), 'flux', gN);

[sol, G] = VEM3D(G, f, bc, k, 'cellProjectors', true);
U = [sol.nodeValues; sol.edgeValues; sol.faceMoments; sol.cellMoments];

Kc = G.cells.centroids;
cells = 1:G.cells.num;
r = .7; c = [1,0,0];
cells = cells(sum(bsxfun(@minus, Kc, c).^2,2) > r^2);
faceNum = mcolon(G.cells.facePos(cells),G.cells.facePos(cells+1)-1);
faces = G.cells.faces(faceNum);

sol = calculateFaceAverages(G,sol);


figure();
plotFaces(G,faces,sol.faceMoments(faces));
colorbar;
view(3);
axis equal;

if k == 2

IF = polygonInt3D(G,1:G.faces.num,gD, 7);
IC = polyhedronInt(G,1:G.cells.num,gD, 7);

u = [gD([G.nodes.coords; G.edges.centroids]); IF./G.faces.areas; IC./G.cells.volumes];
err = abs((U - u));
elseif k == 1
    u = gD(G.nodes.coords);
    err = abs(U-u);
end



l2Err = l2Error3D(G, sol, gD ,k);

fprintf('2-norm error: %d\n', norm(err, 2));
fprintf('L^2-norm error: %d\n\n', sqrt(sum(l2Err.^2)));
figure()
plot(err);