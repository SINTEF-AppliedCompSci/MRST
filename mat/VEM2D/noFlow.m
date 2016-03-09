clc; clear; close all;
addpath('../')

nx = 30; ny = 30;
xMax = 1; yMax = 1;
G = cartGrid([nx, ny], [xMax, yMax]);

G = computeGeometry(G);
G = mrstGridWithFullMappings(G);
% G = sortEdges(G);


bcIn = find(G.faces.centroids(:,1) == 0);
bcNoFlow = find(G.faces.centroids(:,1) == xMax/2 & G.faces.centroids(:,2) < .7*yMax);
bcOut = find(G.faces.centroids(:,1) == xMax & G.faces.centroids(:,2) > .5*yMax);
bcDir = find(G.faces.centroids(:,2) == 0 | G.faces.centroids(:,2) == yMax);

gNIn = @(X) -ones(size(X,1),1);
gNOut = @(X) ones(size(X,1),1);
gD = @(X) zeros(size(X,1),1);

bc1 = struct('bcFunc', {{gNIn, gD, gNOut, gD}}, 'bcFaces', {{bcIn, bcNoFlow, bcOut, bcDir}}, 'bcType', {{'neu', 'dir','neu','dir'}});
bc2 = struct('bcFunc', {{gNIn, gNOut, gD}}, 'bcFaces', {{bcIn, bcOut, bcDir}}, 'bcType', {{'neu','neu','dir'}});


f = @(X) zeros(size(X,1),1);

U1 = VEM2D(G,f,bc1);
U2 = VEM2D(G,f,bc2);

subplot(1,2,1)
plotVEM(G, U1, '')
X = G.nodes.coords;

subplot(1,2,2)
plotVEM(G, U2, '')

% k = 1;
% alpha = 1;
% 
% G = computeVEM2DGeometry(G,f,k, alpha);
% 
% [A,b] = VEM2D_glob_v2(G, bc1, k);
% 
% N = G.nodes.num + G.faces.num*(k-1) + G.cells.num*k*(k-1)/2;
% 
% [bcDof, bBC] = VEM2D_bc(G,bc1,k);
% b(bcDof == 1) = bBC(bcDof == 1);
% SBC = spdiags(ones(N,1),0,N,N);
% A(bcDof == 1,:) = SBC(bcDof == 1,:);
% b(bcDof == 2) = b(bcDof == 2) + bBC(bcDof == 2);
% 
% U = A\b;