clc; clear; close all;

%   TEST 4: Finite difference 3D

nx = 3; ny = 3; nz = 3;
G = cartGrid([nx,ny,nz],[1,ny/nx,nz/nx]);
h = sqrt(3)/nx;

% beta =1.8395265*10e-5;
beta = -97.2;
w1 = beta*1/4; w2 = (1-beta)*9/10+beta*3/4; w3 = 1-w1-w2;

[A_FD, epsx, epsy, epsz] = stencils3D(G,w1,w2,w3);



alpha = (1/20*beta +1/5)*sqrt(3)*h;

f = @(X) zeros(size(X,1),1);
G = computeVEMGeometry(G,f,1);

boundaryEdges = find((G.faces.neighbors(:,1) == 0) + (G.faces.neighbors(:,2) == 0));
gD = @(X) X(:,3);
bc = struct('bcFunc', {{gD}}, 'bcFaces', {{boundaryEdges}}, 'bcType', {{'dir'}});
nN = G.nodes.num;
[bcDof, bBC] = VEM3D_bc(G,bc,1);
SBC = spdiags(ones(nN,1),0,nN,nN);

[sol, A_VEM,b_VEM] = VEM3D(G,f,bc, 1, alpha*ones(G.cells.num,1));
A_FD(bcDof == 1,:) = SBC(bcDof == 1,:);

subplot(1,2,1);
spy(A_FD);
subplot(1,2,2);
spy(abs(A_VEM)>10e-10);

diff = norm(A_VEM-A_FD,'fro')

% fChi = f(G.nodes.coords);
% b_FEM = epsx*epsy*epsz*fChi;
% b_FEM(bcDof == 1) = gD(G.nodes.coords(bcDof == 1, :));
% U_FD = A_FD\b_FEM;
% u = gD(G.nodes.coords);
% norm(u - U_FD,2);
% 
% U_VEM = A_VEM\b_VEM;

