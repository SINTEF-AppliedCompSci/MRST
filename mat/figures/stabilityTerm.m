clc; clear; close all;

addpath('../VEM2D');

%%

nx = 151; ny = 11;
xMax = 2; yMax = 2;
C = [xMax,yMax]/2;

% C = -[.05, .05];

G = cartGrid([nx,ny], [xMax,yMax]);
G.nodes.coords = bsxfun(@minus, G.nodes.coords,C);
G = computeVEM2DGeometry(G);

%%

diff = sum(G.cells.centroids.^2,2);

sourceCell = find(diff == min(diff));
src = addSource([],sourceCell(1), 1);

C = [0,0];
gD = @(X) -1/(2*pi)*log( sqrt( sum(bsxfun(@minus, X, C).^2,2) ) );

bF = find(any(G.faces.neighbors == 0,2));
bc = VEM2D_addBC([], G, bF, 'pressure', gD);

%%

[sol1,G] = VEM2D(G, 0, bc, 1,'projectors', true, 'src', src);

%%

X = G.nodes.coords(G.cells.nodes(G.cells.nodePos(1):G.cells.nodePos(2)-1),:);
hx = abs(max(X(:,1))-min(X(:,1)));
hy = abs(max(X(:,2))-min(X(:,2)));
Q = sqrt(9/(4*hx*hy))*[-1,1,-1,1]';

lambda = 9*1.38*(1/hx^2+1/hy^2);
sigma = 1/(Q'*Q)^2*lambda*ones(G.cells.num,1);

[sol2,G] = VEM2D(G, 0, bc, 1, 'projectors', true, 'sigma', sigma, 'cartGridQ',  true, 'src', src);

%%

azel = [-19,14];
uChi = gD(G.nodes.coords);
dest = '../../tex/thesis/fig/';

fig1 = figure;
plotVEM2D(G, sol1, 1)
set(gcf, 'defaultTextInterpreter', 'LaTex');
set(gca,'XTick', [-1,0,1]);
set(gca,'YTick', [-1,0,1]);
xlabel('$x$'); ylabel('$y$'); zlabel('$u_h$');
view(azel);

savePdf(gcf,strcat(dest,'stabilityBad.pdf'))

fig2 = figure;
plotVEM2D(G, sol2, 1)
set(gcf, 'defaultTextInterpreter', 'LaTex');
set(gca,'XTick', [-1,0,1]);
set(gca,'YTick', [-1,0,1]);
xlabel('$x$'); ylabel('$y$'); zlabel('$u_h$');
view(azel);

savePdf(gcf,strcat(dest, 'stabilityGood.pdf'))

fig3 = figure;
sol = sol1;
sol.nodeValues = uChi;
plotVEM2D(G, sol, 1)
view(azel);

%%
r = .1;
inside = sum(G.nodes.coords.^2,2) < r^2;

err = abs(sol1.nodeValues-uChi);
err(inside) = max(err(~inside));
sol = sol1;
sol.nodeValues = err;
sol = calculateCellAverages(G,sol);
fig4 = figure;
plotCellData(G, sol.cellMoments)
colorbar;

err = abs(sol2.nodeValues-uChi);
err(inside) = max(err(~inside));
sol = sol2;
sol.nodeValues = err;
sol = calculateCellAverages(G,sol);
fig5 = figure;
plotCellData(G, sol.cellMoments)
colorbar;




%%

l2err1 = sqrt(sum(l2Error(G, sol1, gD, 1)));
err1 = norm((sol1.nodeValues-uChi),2);
l2err2 = sqrt(sum(l2Error(G, sol2, gD, 1)));
err2 = norm((sol2.nodeValues-uChi),2);

epsilon = hx/hy;

fprintf('Aspect ratio: %f \n', epsilon);
fprintf('Error with default stability term: %f, \n', err1);
fprintf('L^2-Error with default stability term: %f, \n', l2err1);
fprintf('Error with custom stability term: %f, \n', err2);
fprintf('L^2-Error with default stability term: %f, \n', l2err2);

%%
% 
% X = G.nodes.coords(G.cells.nodes(G.cells.nodePos(1):G.cells.nodePos(2)-1),:);
% hx = abs(max(X(:,1))-min(X(:,1)));
% hy = abs(max(X(:,2))-min(X(:,2)));
% Q = sqrt(9/(4*hx*hy))*[-1,1,-1,1]';
% 
% uChi = gD(G.nodes.coords);
% 
% wLower = 0; wUpper = 10;
% tol = 1e-2;
% while abs(wLower - wUpper) > tol;
% 
% n = 10;
% w = linspace(wLower,wUpper,n);
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
% i = find(min(errVec) == errVec);
% wMin = w(i);
% iLower = max(0,i-1);
% wLower = w(iLower);
% iUpper = min(n, i+1);
% wUpper = w(iUpper);
% 
% end
% 
% 
% lambda = 9*wMin*(1/hx^2+1/hy^2);
% sigma = 1/(Q'*Q)^2*lambda*ones(G.cells.num,1);
% fprintf('Minimum attained at w = %f', wMin);
% fprintf('Minimum attained at sigma = %f', sigma(1));