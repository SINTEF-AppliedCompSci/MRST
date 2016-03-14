clc; clear; close all;

addpath('../'); addpath('../VEM3D/');

%   TEST 4: Finite difference 3D

nx = 2; ny = 2; nz = 2;
G = cartGrid([nx,ny,nz]);
h = sqrt(3)/nx;

% % beta =1.8395265*10e-5;
% beta = -97.2;
% w1 = beta*1/4; w2 = (1-beta)*9/10+beta*3/4; w3 = 1-w1-w2;
% 


w1 = 1/99; w2 = 90/8; w3 = 1-w1-w2;
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

% ai = 0:.1:2; ni = numel(ai);
% aj = 0:.1:6; nj = numel(aj);
% ak = 0:.25:1; nk = numel(ak);
% al = 0:.25:1; nl = numel(al);
% out = [];
% for i = 1:ni
%     for j = 1:nj
%         alpha = [6*ai(i), 6*ai(i), 6*ai(i), 9*aj(j)];
%         alphaMat = repmat(alpha, G.cells.num, 1);
%         [~,A_VEM,b_VEM, ~] = VEM3D(G,f,bc,1,alphaMat);
%         out = [out ; ai(i), ai(i), ai(i) aj(j), norm(A_VEM-A_FD,'fro')];
%         alpha
%     end
% end

alpha = [6*(3/2*w1 + w2), 6*(3/2*w1 + w2), 6*(3/2*w1 + w2), 9*(9/2*w1 + 3/2*w2 + 9/2*w3)];
alpha = repmat(alpha,G.cells.num, 1);
[~,A_VEM,b_VEM, ~] = VEM3D(G,f,bc,1,alpha);
norm(A_VEM-A_FD, 'fro')


%     

% for i = 1:ni
%     alpha = [ai(i), ai(i), ai(i), ai(i), 1,1,1,0,0,0];
%     alphaMat = h*repmat(alpha,G.cells.num,1);
%     [sol, A_VEM,b_VEM] = VEM3D(G,f,bc, 1, alphaMat);
%     out = [out ; alpha, norm(A_VEM-A_FD,'fro')];
%     alpha
% end

%     alpha = 1/8*[24, 24, 24, 24, 1,1,1,0,0,0];
%     alphaMat = h*repmat(alpha,G.cells.num,1);
%     [sol, A_VEM,b_VEM] = VEM3D(G,f,bc, 1, alphaMat);
%     out = [out ; alpha, norm(A_VEM-A_FD,'fro')];


% minimum = out((min(out(:,end)) == out(:,end)),:);
% 
% plot(out(:,end))
% 
% fprintf('Minimal difference %f obtained with alpha = (%f, %f, %f, %f)', minimum(5), minimum(1), minimum(2), minimum(3), minimum(4))
