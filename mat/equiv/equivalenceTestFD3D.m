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


w1 = 1/4; w2 = 3/4; w3 = 0;
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

ai = 0:.25:1; ni = numel(ai);
aj = 0:.25:1; nj = numel(aj);
ak = 0:.25:1; nk = numel(ak);
al = 0:.25:1; nl = numel(al);
am = 0:.25:1; nm = numel(am);
an = 0:.25:1; nn = numel(an);
ao = 0:.25:1; no = numel(ao);

out = [];

for i = 1:ni
    for j = 1:nj
        for k = 1:nk
            for l = 1:nl
                for m = 1:nm
                    for n = 1:nn
                        for o = 1:no
                            alpha = [ai(i), ai(i), ai(i), ai(i), aj(j), ak(k), al(l), am(m), an(n), ao(o)];
                            alphaMat = repmat(alpha,G.cells.num,1);
                            [~, A_VEM,~] = VEM3D(G,f,bc, 1, alphaMat);
                            out = [out ; alpha, norm(A_VEM-A_FD,'fro')];
                        end
                    end
                end
            end
        end
    end
    alpha
end
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


minimum = out((min(out(:,end)) == out(:,end)),:);

plot(out(:,end))

fprintf('Minimal difference %f obtained with alpha = (%f, %f, %f, %f)', minimum(5), minimum(1), minimum(2), minimum(3), minimum(4))
