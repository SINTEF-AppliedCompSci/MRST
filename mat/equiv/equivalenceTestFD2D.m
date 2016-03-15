clc; clear; close all;

addpath('../'); addpath('../VEM2D/');

%   TEST 3: Finite difference
nx = 5; ny = 5;
xMax = 5; yMax = 5;
G = cartGrid([nx,ny],[xMax, yMax]);
f = @(X) zeros(size(X,1),1);
gD = @(X) X(:,1);

boundaryEdges = find((G.faces.neighbors(:,1) == 0) + (G.faces.neighbors(:,2) == 0));
bc = struct('bcFunc', {{gD}}, 'bcFaces', {{boundaryEdges}}, 'bcType', {{'dir'}});

w = 99;
hx = xMax/(2*nx); hy = yMax/(2*ny);
epsilon = hy/hx;
alpha = 3*(1/hx^2 + 1/hy^2)*w*2;
k = 1;
G = computeVEM2DGeometry(G,f,k,alpha);

[A,b] = VEM2D_glob_v2(G, bc, k);

nN = G.nodes.num;


main = -4*ones(nN,1);
sub = ones(nN,1);

A_FD_V = -epsilon*spdiags([sub,.5*main,sub], [-1,0,1], nN, nN) ...
         -1/epsilon*spdiags([sub,.5*main,sub], [-(nx+1),0,nx+1], nN, nN);
% A_FD_V = -spdiags([sub,sub,main,sub,sub], [-sqrt(nN), -1,0,1, sqrt(nN)], nN, nN);
A_FD_D = -.5*spdiags([sub,sub,main,sub,sub], [-(nx+2), -nx, 0, nx, nx+2], nN, nN);

wD =1-w; wV = w;
A_FD = wD*A_FD_D + wV*A_FD_V;

[bcDof, bBC] = VEM2D_bc(G,bc,k);
SBC = spdiags(ones(nN,1),0,nN,nN);
A(bcDof == 1,:) = SBC(bcDof == 1,:);
A_FD(bcDof == 1,:) = SBC(bcDof == 1,:);

diff = norm(A-A_FD,'fro')
% full(A)
% full(A_FD)