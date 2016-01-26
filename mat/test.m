clc; clear all; close all;

run('../../matlab/project-mechanics-fractures/mystartup.m')

n = 15;

% G = cartGrid([2,2,30],[1,1,30]);
G = cartGrid([n,n,n],[1,1,1]);
% G.nodes.coords(1:2,:) = G.nodes.coords(1:2,:) -0.5;
G = computeGeometry(G);
G = mrstGridWithFullMappings(G);

f = @(X) zeros(size(X,1),1);

% f = @(X) X(:,3);


G = computeVEMGeometry(G,f);

boundaryFaces = (1:G.faces.num)';

boundaryFaces = boundaryFaces( G.faces.centroids(:,1) == 0 | ...
                               G.faces.centroids(:,1) == 1 | ...
                               G.faces.centroids(:,2) == 0 | ...
                               G.faces.centroids(:,2) == 1 | ...
                               G.faces.centroids(:,3) == 0 | ...
                               G.faces.centroids(:,3) == 1 );                          

% boundaryFaces = boundaryFaces(G.faces.centroids(:,3) == 0);
                         
C = -[.2,.2,.2];
gD = @(X) -1./(2*pi*sqrt(sum((X-repmat(C,size(X,1),1)).^2,2)));
            
bc = struct('bcFunc', {{gD}}, 'bcFaces', {{boundaryFaces}}, 'bcType', {{'dir'}});

[A,b] = VEM3D_glob(G,f,bc);

U = A\b;

nodeValues = full(U(1:G.nodes.num));
edgeMidValues = full(U(G.nodes.num + 1:G.nodes.num + G.edges.num));
faceAvg = full(U(G.nodes.num + G.edges.num + 1: G.nodes.num + G.edges.num + G.faces.num));
cellAvg = full(U(G.nodes.num + G.edges.num + G.faces.num + 1:end));

figure();
plotFaces(G, 1:G.faces.num, faceAvg);
colorbar;
view(3);

err = U - gD([G.nodes.coords; G.edges.centroids; G.faces.centroids; G.cells.centroids])
% plot(nodeValues);

%   Implement: efficient rule for faceIntegrals.
%              change bc to give avg values.