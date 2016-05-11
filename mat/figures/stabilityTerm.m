clc; clear; close all;

addpath('../VEM2D');

%%

nx = 151; ny = 11;
xMax = 2; yMax = 2;
C = [xMax,yMax]/2;

G = cartGrid([nx,ny], [xMax,yMax]);
G.nodes.coords = bsxfun(@minus, G.nodes.coords,C);
G = computeVEM2DGeometry(G);

%%

diff = sum(G.cells.centroids.^2,2);

sourceCell = find(diff == min(diff));
src = addSource([],sourceCell(1), 1);
gD = @(X) -1/(2*pi)*log( sqrt( sum(X.^2,2) ) );

bF = find(any(G.faces.neighbors == 0,2));
bc = VEM2D_addBC([], G, bF, 'pressure', gD);

%%

[sol1,G] = VEM2D(G, 0, bc, 1,'projectors', true, 'src', src);

%%

X = G.nodes.coords(G.cells.nodes(G.cells.nodePos(1):G.cells.nodePos(2)-1),:);
hx = abs(max(X(:,1))-min(X(:,1)));
hy = abs(max(X(:,2))-min(X(:,2)));
Q = sqrt(9/(4*hx*hy))*[-1,1,-1,1]';

lambda = 9*1.39*(1/hx^2+1/hy^2);
sigma = 1/(Q'*Q)^2*lambda*ones(G.cells.num,1);

[sol2,G] = VEM2D(G, 0, bc, 1, 'projectors', true, 'sigma', sigma, 'cartGridQ',  true, 'src', src);

%%

azel = [-19,14];
uChi = gD(G.nodes.coords);
dest = '../../tex/thesis/fig/';

fig1 = figure;
plotVEM2D(G, sol1, 1)
set(gcf, 'defaultTextInterpreter', 'LaTex');
xlabel('$x$'); ylabel('$y$'); zlabel('$u_h$');
view(azel);

savePdf(gcf,strcat(dest,'stabilityBad.pdf'))

fig2 = figure;
plotVEM2D(G, sol2, 1)
set(gcf, 'defaultTextInterpreter', 'LaTex');
xlabel('$x$'); ylabel('$y$'); zlabel('$u_h$');
view(azel);

savePdf(gcf,strcat(dest, 'stabilityGood.pdf'))

err1 = sqrt(sum(l2Error(G, sol1, gD, 1)));
err2 = sqrt(sum(l2Error(G, sol2, gD, 1)));
err2 = norm((sol2.nodeValues-uChi)./uChi,2);

epsilon = hx/hy;

fprintf('Aspect ratio: %f \n', epsilon);
fprintf('Error without default stability term: %f, \n', err1);
fprintf('Error with default stability term: %f, \n', err2);

%%
% 
% n = 10;
% w = linspace(1.391,1.395,n);
% errVec = zeros(n,1);
% 
% for i = 1:n
%     
%     lambda = 9*w(i)*(1/hx^2+1/hy^2);
%     sigma = 1/(Q'*Q)^2*lambda*ones(G.cells.num,1);
% 
%     [sol2,G] = VEM2D(G, 0, bc, 1, 'projectors', true, 'sigma', sigma, 'cartGridQ',  true, 'src', src);
%     
%     errVec(i) = norm(sol2.nodeValues-uChi,2);
%     
% end
% 
% i = find(min(errVec) == errVec)
% w(i)

