%--------------------------------------------------------------------------
%   Minimal working example of VEM in three dimensions.
%-----------------------------------------------------------------ØSK-2016-

clc; clear; close all;

addpath('../')              %  Extra grid mappings.
            
                            %   Voronoi grid generation. Replace path with
                            %   your own local version.
addpath('/home/strene/Documents/master/coop/pebiGridding/voronoi3D');

ex = 1;

switch ex
    case 1
        k = 2;
        gD = @(X) X(:,1) + X(:,3)/77 - X(:,2);
        f = 0;
        gridType = 'cart';
        n = 100;
        nx = 4; ny = 4; nz = 4;
        xMax = 1; yMax = 1; zMax = 1;
end

if strcmp(gridType, 'cart')
    G = cartGrid([nx,ny,nz], [xMax, yMax, zMax]);
else
    G = voronoiCube(n,[1,1,1]);
end

G = computeVEM3DGeometry(G);
G = VEM3D_faceProjectors(G,k);

bFaces = find( any(G.faces.neighbors == 0,2)) ;
bc = VEM3D_addBC([], bFaces, 'pressure', gD);

[sol, G] = VEM3D(G,f,bc,k, 'cellProjectors', true);

Kc = G.cells.centroids;
cells = 1:G.cells.num;
r = .7;
cells = cells(Kc(:,1).^2 + Kc(:,2).^2 + Kc(:,3).^2 > r^2);
faceNum = mcolon(G.cells.facePos(cells),G.cells.facePos(cells+1)-1);
faces = G.cells.faces(faceNum);

if k == 2
    figure();
    plotFaces(G,faces,sol.faceMoments(faces));
    colorbar;
    view(3);
    axis equal;
end

l2Err = l2Error3D(G, sol, gD, k);
if k == 1
    u = gD(G.nodes.coords);
    err = sol.nodeValues - u;
    n2Err = norm(err);
elseif k == 2
    u = [gD([G.nodes.coords; G.edges.centroids]); ...
         polygonInt3D(G,1:G.faces.num,gD, 7)./G.faces.areas; ...
         polyhedronInt(G, 1:G.cells.num,gD,7)./G.cells.volumes];
    err = [sol.nodeValues; sol.edgeValues; sol.faceMoments; sol.cellMoments] - u;
    n2Err = norm(err);
end

fprintf('L2-error: \t %d \n\n', sqrt(sum(l2Err.^2)));
fprintf('2-norm of error: \t %d \n\n', n2Err);