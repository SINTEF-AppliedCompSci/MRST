clc; clear all; close all;

run('../../matlab/project-mechanics-fractures/mystartup.m')

n = 15;
nx = n; ny = n;
G = unitSquare(nx, ny);
% G = cartGrid([nx,ny],[1,1]);
G = sortEdges(G);
G = mrstGridWithFullMappings(G);
G = computeGeometry(G);

f = @(X) 4*pi^2*sin(X(:,1)*pi).*cos(X(:,2)*pi);
f = @(X) X(:,1).^2 + 100*X(:,2).^2 + 39/7*X(:,1) + X(:,1).*X(:,2);
gD = @(X) 2*sin(X(:,1)*pi).*cos(X(:,2)*pi) - log(1./((X(:,1)+0.1).^2 + (X(:,2)+0.1).^2));
gN = @(X) -2*pi*cos(X(:,1)*pi).*cos(X(:,2)*pi) - 2*(X(:,1)+0.1)./((X(:,1)+0.1).^2 + (X(:,2)+0.1).^2);

tol = 1e-6;
boundaryEdges = find((G.faces.neighbors(:,1) == 0) + (G.faces.neighbors(:,2) == 0));
bNeu = abs(G.faces.centroids(boundaryEdges,1)) < tol;
bc = struct('bcFunc', {{gN, gD}}, 'bcFaces', {{boundaryEdges(bNeu), boundaryEdges(~bNeu)}}, 'bcType', {{'neu', 'dir'}});

% gD = @(X) 1/2.*(X(:,2).^2 - X(:,1).^2);
% gN = @(X) X(:,2);
% f = @(X) zeros(size(X,1),1);
% boundaryEdges = find((G.faces.neighbors(:,1) == 0) + (G.faces.neighbors(:,2) == 0));
% bNeu = zeros(numel(boundaryEdges),1);
% bDir = zeros(numel(boundaryEdges),1);
% for e = 1:numel(boundaryEdges)
%     X = G.faces.centroids(boundaryEdges(e),:);
%     if X(:,1) < 0 + eps 
%         bDir(e) = 1;
%     end
%     if X(:,1) > 1 - eps
%         bDir(e) = 1;
%     end
%     if X(:,2) < 0 + eps 
%         bNeu(e) = 1;
%     end
%     if X(:,2) > 1-eps
%         bNeu(e) = 1;
%     end
%     
% end
% 
% bc = struct('bcFunc', {{gN,gD}}, 'bcFaces', {{boundaryEdges(bNeu == 1), boundaryEdges(bDir == 1)}}, 'bcType', {{'neu', 'dir'}});


U = VEM2D(G,f,bc);

X = [G.nodes.coords ; G.faces.centroids];
Uexact = gD(X);

errVec = U(1:end-G.cells.num) - Uexact;
%errVec = U - Uexact;
err = norm(errVec, inf);
fprintf('Error: %e', err);

% fig1 = figure;
% plotGridWithDofs(G,bc);

fig2 = figure;
plotVEM(G, U, '');

fig5 = figure;

% fig3 = figure;
% plotVEM(G, U, 'dof');

% fig4 = figure;
% plot(errVec);

%   Make function for computing baricenters of all cells
%   Make function for detecting boundary dofs
%   Change VEM2D to accept force term f