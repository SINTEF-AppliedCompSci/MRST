clc; clear; close all;

run('../../matlab/project-mechanics-fractures/mystartup.m')

f = @(X) ones(size(X,1),1);
gD = @(X) -(X(:,1).^2 + X(:,2).^2 + X(:,3).^2)/6;

nVec = [2, 4, 6, 8];
hVec = zeros(1,4);

gridLim = [1,1,1];

for i = 1:4
    
    n = nVec(i);
    G = cartGrid([n,n,n],gridLim);
%     G = unitCubeTetrahedrons(n);

    G = computeGeometry(G);
    G = mrstGridWithFullMappings(G);
    G = computeVEMGeometry(G,f);

    boundaryFaces = (1:G.faces.num)';
    boundaryFaces = boundaryFaces( G.faces.centroids(:,1) == 0          | ...
                                   G.faces.centroids(:,1) == gridLim(1) | ...
                                   G.faces.centroids(:,2) == 0          | ...
                                   G.faces.centroids(:,2) == gridLim(2) | ...
                                   G.faces.centroids(:,3) == 0          | ...
                                   G.faces.centroids(:,3) == gridLim(3) );                          

    bc = struct('bcFunc', {{gD}}, 'bcFaces', {{boundaryFaces}}, 'bcType', {{'dir'}});

    [A,b] = VEM3D_glob(G,f,bc);

    U = A\b;

    IF = polygonInt3D(G,1:G.faces.num,gD);
    IC = polyhedronInt(G,1:G.cells.num,gD);

    u = [gD([G.nodes.coords; G.edges.centroids]); IF./G.faces.areas; IC./G.cells.volumes];
    
    h = 0;
    nN = G.cells.num;
    for j = 1:nN
        nodeNum = G.cells.nodePos(j):G.cells.nodePos(j+1)-1;
        nodes = G.cells.nodes(nodeNum);
        X = G.nodes.coords(nodes,:);
        h = h + cellDiameter(X);
    end
    h = h/nN;
    hVec(i) = h;
    errVec(i) = h^(3/2)*norm(U - u, inf);
end


loglog(hVec, errVec);

slope = polyfit(log(hVec), log(errVec), 1);

    