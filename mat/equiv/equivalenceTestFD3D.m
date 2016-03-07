clc; clear; close all;

addpath('../'); addpath('../VEM3D/');

%   TEST 4: Finite difference 3D

nx = 2; ny = 2; nz = 2;
G = cartGrid([nx,ny,nz],[1,ny/nx,nz/nx]);
h = sqrt(3)/nx;

% % beta =1.8395265*10e-5;
% beta = -97.2;
% w1 = beta*1/4; w2 = (1-beta)*9/10+beta*3/4; w3 = 1-w1-w2;
% 


w1 = 1; w2 = 0; w3 = 0;
[A_FD, epsx, epsy, epsz] = stencils3D(G,w1,w2,w3);

f = @(X) zeros(size(X,1),1);
G = computeVEMGeometry(G,f,1);

boundaryEdges = find((G.faces.neighbors(:,1) == 0) + (G.faces.neighbors(:,2) == 0));
gD = @(X) X(:,3);
bc = struct('bcFunc', {{gD}}, 'bcFaces', {{boundaryEdges}}, 'bcType', {{'dir'}});
nN = G.nodes.num;
[bcDof, bBC] = VEM3D_bc(G,bc,1);
SBC = spdiags(ones(nN,1),0,nN,nN);
h = G.cells.diameters(1);

A_FD(bcDof == 1,:) = SBC(bcDof == 1,:);

v = .62:.005:.78;
nv = size(v,2);
alphaMat = [rldecode(v,nv^3*ones(1,nv),2); ...
	 repmat(rldecode(v,nv^2*ones(1,nv),2),1,nv); ...
	 repmat(rldecode(v,nv*ones(1,nv),2),1,nv^2); ...
	 repmat(v,1,nv^3)];
nAlpha = size(alphaMat,2);

for i = 1:nAlpha
	  alpha = h*repmat(alphaMat(:,i)',G.cells.num,1);
            [sol, A_VEM,b_VEM] = VEM3D(G,f,bc, 1, alpha);
            diff(i) = norm(A_VEM-A_FD,'fro');
end
    

alpha = alphaMat(:,find(min(diff) == diff));
diffMin = min(diff);

fprintf('Minimal difference %f obtained with alpha = (%f, %f, %f, %f)', diffMin, alpha(1), alpha(2), alpha(3), alpha(4))


% subplot(1,2,1);
% spy(A_FD);
% subplot(1,2,2);
% spy(abs(A_VEM)>10e-10);



% fChi = f(G.nodes.coords);
% b_FEM = epsx*epsy*epsz*fChi;
% b_FEM(bcDof == 1) = gD(G.nodes.coords(bcDof == 1, :));
% U_FD = A_FD\b_FEM;
% u = gD(G.nodes.coords);
% norm(u - U_FD,2);
% 
% U_VEM = A_VEM\b_VEM;

