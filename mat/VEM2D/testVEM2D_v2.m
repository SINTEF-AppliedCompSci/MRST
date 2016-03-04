clc; clear all; close all;

addpath('../');


G = cartGrid([10,10], [1,1]);
% G = unitSquareTri([4,4],[1,1]);

f = @(X) 0*ones(size(X,1),1);

gD = @(X) 100*X(:,1).^2 - 100*X(:,2).^2;

boundaryEdges = find((G.faces.neighbors(:,1) == 0) + (G.faces.neighbors(:,2) == 0));
bc = struct('bcFunc', {{gD}}, 'bcFaces', {{boundaryEdges}}, 'bcType', {{'dir'}});

k = 1;
alpha = 1;

G = computeVEM2DGeometry(G,f,k, alpha);

[A,b] = VEM2D_glob_v2(G, bc, k);

N = G.nodes.num + G.faces.num*(k-1) + G.cells.num*k*(k-1)/2;

[bcDof, bBC] = VEM2D_bc(G,bc,k);
b(bcDof == 1) = bBC(bcDof == 1);
SBC = spdiags(ones(N,1),0,N,N);
A(bcDof == 1,:) = SBC(bcDof == 1,:);
b(bcDof == 2) = b(bcDof == 2) + bBC(bcDof == 2);

U = A\b;

if k == 1
    u = gD(G.nodes.coords);
elseif k == 2
    gDint = polygonInt_v2(G, 1:G.cells.num, gD, 7);
    u = [gD([G.nodes.coords; G.faces.centroids]); gDint./G.cells.volumes];
end

% plotVEM(G,U,'dof')    

err = U-u;

norm(err,inf)


