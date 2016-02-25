clc; clear all; close all;
addpath('../');
addpath('../VEM2D/');

G = cartGrid([1,1], [2,2]);

X = G.nodes.coords;
X = bsxfun(@minus, X, [1,1]);
G.nodes.coords = X;

f = @(X) X(:,1);


% %   TEST 1: unit square
% 
% 
% k = 1;
% 
% G = computeVEM2DGeometry(G,f,k,2/3);
% 
% load('FEM2D_1st_sq');
% AK_VEM = G.cells.AK{1};
% 
% err = norm(AK_FEM - AK_VEM, 'fro');
% 
% fprintf('Result from unit square comparison: \t  %d \n\n', err);

%   TEST 2: unit suqare 2nd

k = 2;
G = computeVEM2DGeometry(G,f,k,0);
AK_VEM = G.cells.AK{1};
SK = G.cells.SK{1};
PN = G.cells.PN{1};



load('FEM2D_2nd_sq');

load('FEM2VEM.mat');

Fi= inv(F);
AK_FEM = F*AK_FEM*F';

% AK_FEM-AK_VEM-SK

alpha = (AK_FEM(1,1) - AK_VEM(1,1))/SK(1,1);


%alpha = 648/95



AK_FEM - AK_VEM - alpha*SK
% abs(AK_FEM - AK_VEM -alpha*SK) < 10e-4



err = norm(AK_FEM - AK_VEM - alpha*SK, 'fro');

fprintf('Result from unit square comparison: \t  %d \n\n', err);

% %   TEST 3: Finite difference
% G = cartGrid([3,3],[1,1]);
% f = @(X) zeros(size(X,1),1);
% gD = @(X) X(:,1);
% 
% tol = 1e-6;
% boundaryEdges = find((G.faces.neighbors(:,1) == 0) + (G.faces.neighbors(:,2) == 0));
% bc = struct('bcFunc', {{gD}}, 'bcFaces', {{boundaryEdges}}, 'bcType', {{'dir'}});
% 
% w = 99/100;
% alpha = 2*w;
% k = 1;
% G = computeVEM2DGeometry(G,f,k,alpha);
% h = G.nodes.coords(2,1)-G.nodes.coords(1,1);
% 
% [A,b] = VEM2D_glob_v2(G, bc, k);
% 
% nN = G.nodes.num;
% 
% 
% main = -4*ones(nN,1);
% sub = ones(nN,1);
% A_FD_V = -spdiags([sub,sub,main,sub,sub], [-sqrt(nN), -1,0,1, sqrt(nN)], nN, nN);
% A_FD_D = -.5*spdiags([sub,sub,main,sub,sub], [-(sqrt(nN)+1), -(sqrt(nN)-1), 0, sqrt(nN)-1, sqrt(nN) + 1], nN, nN);
% 
% wD =1-w; wV = w;
% A_FD = wD*A_FD_D + wV*A_FD_V;
% 
% [bcDof, bBC] = VEM2D_bc(G,bc,k);
% SBC = spdiags(ones(nN,1),0,nN,nN);
% A(bcDof == 1,:) = SBC(bcDof == 1,:);
% A_FD(bcDof == 1,:) = SBC(bcDof == 1,:);
% 
% norm(A-A_FD,'fro')
% full(A)
% full(A_FD)
