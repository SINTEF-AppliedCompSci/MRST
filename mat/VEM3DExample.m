clc; clear; close all;

run('../../matlab/project-mechanics-fractures/mystartup.m')

alpha = 2;
f = @(X) sin(X(:,1)).*cos(alpha*pi*X(:,2)).*X(:,3)*(1+pi^2*alpha^2);
gD = @(X) sin(X(:,1)).*cos(alpha*pi*X(:,2)).*X(:,3);
gridLim = [1,1,1];


nVec = [4, 7, 9, 12];
hVec = zeros(1,4);
errVec = zeros(0,4);



for i = 1:4
    
    n = nVec(i);
%     G = cartGrid([n,n,n],gridLim);
%     G = unitCubeTetrahedrons(n);
    G = voroniCube(n^3,gridLim);

    G = computeGeometry(G);
    G = mrstGridWithFullMappings(G);
    G = computeVEMGeometry(G,f);

    boundaryFaces = (1:G.faces.num)';
    boundaryFaces = boundaryFaces( G.faces.neighbors(:,1) == 0 | ...
                                   G.faces.neighbors(:,2) == 0 );                         

    bc = struct('bcFunc', {{gD}}, 'bcFaces', {{boundaryFaces}}, 'bcType', {{'dir'}});

    [A,b] = VEM3D_glob(G,f,bc);

    U = A\b;

    IF = polygonInt3D(G,1:G.faces.num,gD);
    IC = polyhedronInt(G,1:G.cells.num,gD);

    u = [gD([G.nodes.coords; G.edges.centroids]); IF./G.faces.areas; IC./G.cells.volumes];

    nK = G.cells.num;
    h = sum(G.cells.diameters)/nK;
    hVec(i) = h;
    errVec(i) = h^(3/2)*norm(U - u, 2);
    
    e = abs(U-u);
    
    err = l2Norm(G,e);
    
    errVec(i) = err;
    
end

slope = polyfit(log(hVec), log(errVec), 1);

loglog(hVec, errVec);

slope(1)



    