clc; clear; close all;

addpath('../VEM2D');

%%

nx = 101; ny = 11;
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

%%  DEFAULT STABILITY TERM

[sol1,G] = VEM2D(G, 0, bc, 1,'projectors', true, 'src', src);

%%  CUSTOM STABILITYTERMS

X = G.nodes.coords(G.cells.nodes(G.cells.nodePos(1):G.cells.nodePos(2)-1),:);
hx = abs(max(X(:,1))-min(X(:,1)));
hy = abs(max(X(:,2))-min(X(:,2)));
Q = sqrt(9/(4*hx*hy))*[-1,1,-1,1]';

%   Optimal l2-error

w = 0.0615;

lambda = 3*w*(1/hx^2+1/hy^2);
sigma = 1/(Q'*Q)^2*lambda*ones(G.cells.num,1);

[sol2,G] = VEM2D(G, 0, bc, 1, 'projectors', true, 'sigma', sigma, 'cartGridQ',  true, 'src', src);

%   Optimal 2-norm error

w = 4.1409;

lambda = 3*w*(1/hx^2+1/hy^2);
sigma = 1/(Q'*Q)^2*lambda*ones(G.cells.num,1);

[sol3,G] = VEM2D(G, 0, bc, 1, 'projectors', true, 'sigma', sigma, 'cartGridQ',  true, 'src', src);

%   FEM

w = 1;
lambda = 3*w*(1/hx^2+1/hy^2);
sigma = 1/(Q'*Q)^2*lambda*ones(G.cells.num,1);

[sol4,G] = VEM2D(G, 0, bc, 1, 'projectors', true, 'sigma', sigma, 'cartGridQ',  true, 'src', src);


%%

azel = [-19,14];
uChi = gD(G.nodes.coords);
dest = '../../tex/thesis/fig/';
% 
% fig1 = figure;
% plotVEM2D(G, sol1, 1)
% set(gcf, 'defaultTextInterpreter', 'LaTex');
% set(gca,'XTick', [-1,0,1]);
% set(gca,'YTick', [-1,0,1]);
% xlabel('$x$'); ylabel('$y$'); zlabel('$u_h$');
% view(azel);
% 
% savePdf(gcf,strcat(dest,'stabilityBad.pdf'))
% 
% fig2 = figure;
% plotVEM2D(G, sol2, 1)
% set(gcf, 'defaultTextInterpreter', 'LaTex');
% set(gca,'XTick', [-1,0,1]);
% set(gca,'YTick', [-1,0,1]);
% xlabel('$x$'); ylabel('$y$'); zlabel('$u_h$');
% view(azel);
% 
% savePdf(gcf,strcat(dest, 'stabilityL2Err.pdf'))
% 
% fig3 = figure;
% plotVEM2D(G, sol3, 1)
% set(gcf, 'defaultTextInterpreter', 'LaTex');
% set(gca,'XTick', [-1,0,1]);
% set(gca,'YTick', [-1,0,1]);
% xlabel('$x$'); ylabel('$y$'); zlabel('$u_h$');
% view(azel);
% 
% savePdf(gcf,strcat(dest, 'stability2Err.pdf'))
% 
% fig4 = figure;
% plotVEM2D(G, sol4, 1)
% set(gcf, 'defaultTextInterpreter', 'LaTex');
% set(gca,'XTick', [-1,0,1]);
% set(gca,'YTick', [-1,0,1]);
% xlabel('$x$'); ylabel('$y$'); zlabel('$u_h$');
% view(azel);
% 
% savePdf(gcf,strcat(dest, 'stabilityFEM.pdf'))

%%

%%
X = G.nodes.coords(G.cells.nodes(G.cells.nodePos(1):G.cells.nodePos(2)-1),:);
hx = abs(max(X(:,1))-min(X(:,1)));
hy = abs(max(X(:,2))-min(X(:,2)));
Q = sqrt(9/(4*hx*hy))*[-1,1,-1,1]';

r= hy;
trustCells = find(sum(G.cells.centroids.^2,2) > r^2);
trustNodes = unique(G.cells.nodes(mcolon(G.cells.nodePos(trustCells),G.cells.nodePos(trustCells+1)-1)));

%%

gDl2 = sqrt(sum(polygonInt(G, trustCells, @(X) gD(X).*gD(X), 7)));
gD2 = norm(uChi(trustNodes),2);

l2err1 = sqrt(sum(l2Error(G, sol1, gD, 1)))/gDl2;
err1 = norm((sol1.nodeValues-uChi),2)/gD2;

l2err2 = sqrt(sum(l2Error(G, sol2, gD, 1)))/gDl2;
err2 = norm((sol2.nodeValues-uChi),2)/gD2;

l2err3 = sqrt(sum(l2Error(G, sol3, gD, 1)))/gDl2;
err3 = norm((sol3.nodeValues-uChi),2)/gD2;

l2err4 = sqrt(sum(l2Error(G, sol4, gD, 1)))/gDl2;
err4 = norm((sol4.nodeValues-uChi),2)/gD2;

epsilon = hx/hy;

fprintf('Aspect ratio: %f \n', epsilon);
fprintf('Error with default stability term: %f, %f, %f \n', err1, l2err1, err1+l2err1);
fprintf('Error with l2norm stability term: %f, %f, %f \n', err2, l2err2, err2+l2err2);
fprintf('Error with 2norm stability term: %f, %f, %f \n', err3, l2err3, err3+l2err3);
fprintf('Error with FEM stability term: %f, %f, %f \n', err4, l2err4, err4+l2err4);

%%
X = G.nodes.coords(G.cells.nodes(G.cells.nodePos(1):G.cells.nodePos(2)-1),:);
hx = abs(max(X(:,1))-min(X(:,1)));
hy = abs(max(X(:,2))-min(X(:,2)));
Q = sqrt(9/(4*hx*hy))*[-1,1,-1,1]';

uChi = gD(G.nodes.coords);

wLower = 0; wUpper = 6;
tol = 1e-4;
while abs(wLower - wUpper) > tol;

n = 10;
w = linspace(wLower,wUpper,n);
errVec = zeros(n,1);

for i = 1:n
    
    lambda = 3*w(i)*(1/hx^2+1/hy^2);
    sigma = 1/(Q'*Q)^2*lambda*ones(G.cells.num,1);

    [sol,G] = VEM2D(G, 0, bc, 1, 'projectors', true, 'sigma', sigma, 'cartGridQ',  true, 'src', src);
    
%    err = sol.nodeValues - uChi;    
%     errVec(i) = norm(err(trustNodes),2);
    l2Err = l2Error(G, sol, gD, 1);
    errVec(i) = sqrt(sum(l2Err(trustCells)));
    
end

i = find(min(errVec) == errVec);
i = i(1);
wMin = w(i);
iLower = max(1,i-1);
wLower = w(iLower);
iUpper = min(n, i+1);
wUpper = w(iUpper);

end

lambda = 3*wMin*(1/hx^2+1/hy^2);
sigma = 1/(Q'*Q)^2*lambda*ones(G.cells.num,1);
fprintf('Minimum attained at w = %f', wMin);
fprintf('Minimum attained at sigma = %f', sigma(1));

%%
% 
% X = G.nodes.coords(G.cells.nodes(G.cells.nodePos(1):G.cells.nodePos(2)-1),:);
% hx = abs(max(X(:,1))-min(X(:,1)));
% hy = abs(max(X(:,2))-min(X(:,2)));
% Q = sqrt(9/(4*hx*hy))*[-1,1,-1,1]';
% 
% uChi = gD(G.nodes.coords);
% 
% wLower = 0; wUpper = 3;
% n = 20;
% w = linspace(wLower,wUpper,n);
% errVec = zeros(n,3);
% 
% for i = 1:n
%     
%     fprintf('Simulation %d of %d\n\n', i, n)
%     
%     lambda = 9*w(i)*(1/hx^2+1/hy^2);
%     sigma = 1/(Q'*Q)^2*lambda*ones(G.cells.num,1);
% 
%     [sol2,G] = VEM2D(G, 0, bc, 1, 'projectors', true, 'sigma', sigma, 'cartGridQ',  true, 'src', src);
%     
%     errVec(i,1) = sigma(1);
%     errVec(i,2) = norm(sol2.nodeValues-uChi,2)/gD2;
%     errVec(i,3) = sqrt(sum(l2Error(G, sol2, gD, 1)))/gDl2;
%     
% end
% 
% % %%
% % 
% % plot(errVec(:,1), errVec(:,2), errVec(:,1), errVec(:,3));
% 
% 

%%

% gDl2 = sqrt(sum(polygonInt(G, 1:G.cells.num, @(X) gD(X).*gD(X), 7)));
% 
% X = G.nodes.coords(G.cells.nodes(G.cells.nodePos(1):G.cells.nodePos(2)-1),:);
% hx = abs(max(X(:,1))-min(X(:,1)));
% hy = abs(max(X(:,2))-min(X(:,2)));
% Q = sqrt(9/(4*hx*hy))*[-1,1,-1,1]';
% 
% uChi = gD(G.nodes.coords);
% 
% wLower = 0; wUpper = 10;
% n = 20;
% w = linspace(wLower,wUpper,n);
% errVec = zeros(n,1);
% 
% for i = 1:n
%     
%     lambda = 3*w(i)*(1/hx^2+1/hy^2);
%     sigma = 1/(Q'*Q)^2*lambda*ones(G.cells.num,1);
% 
%     [sol,G] = VEM2D(G, 0, bc, 1, 'projectors', true, 'sigma', sigma, 'cartGridQ',  true, 'src', src);
%     
% %     errVec(i) = norm(sol.nodeValues-uChi,2);
%     l2Err = l2Error(G, sol, gD, 1);
%     errVec(i) = sqrt(sum(l2Err(trustCells)));
%     
% end
% 
% %%
% 
% plot(w, errVec/gDl2);
% 
% set(gcf, 'defaulttextinterpreter', 'latex');
% xlabel('$w$'); ylabel('$\|e_h\|_{0,K}/\|u\|_{0,K}$')
% 
% w = 0;
% h = 0;
% ps = get(gcf, 'Position');
% ratio = 1;
% paperWidth = 10;
% paperHeight = paperWidth*ratio;
% set(gcf, 'paperunits', 'centimeters');
% set(gcf, 'papersize', [paperWidth-w paperHeight-h]);
% set(gcf, 'PaperPosition', [-w    -h   paperWidth+w paperHeight+h]);
% print(gcf, '-dpdf', '../../tex/thesis/fig/wStabilityL2.pdf');

%%

X = G.nodes.coords(G.cells.nodes(G.cells.nodePos(1):G.cells.nodePos(2)-1),:);
hx = abs(max(X(:,1))-min(X(:,1)));
hy = abs(max(X(:,2))-min(X(:,2)));
Q = sqrt(9/(4*hx*hy))*[-1,1,-1,1]';

uChi = gD(G.nodes.coords);

wLower = 0; wUpper = 10;
n = 20;
w = linspace(wLower,wUpper,n);
errVec = zeros(n,2);

for i = 1:n
    
    lambda = 3*w(i)*(1/hx^2+1/hy^2);
    sigma = 1/(Q'*Q)^2*lambda*ones(G.cells.num,1);

    [sol,G] = VEM2D(G, 0, bc, 1, 'projectors', true, 'sigma', sigma, 'cartGridQ',  true, 'src', src);
    err = sol.nodeValues - uChi;
    errVec(i,1) = norm(err(trustNodes),2);
    l2Err = l2Error(G, sol, gD, 1);
    errVec(i,2) = sqrt(sum(l2Err(trustCells)));
    
end

%%

plot(w, errVec(:,1)/gD2, w, errVec(:,2)/gDl2);

set(gcf, 'defaulttextinterpreter', 'latex');
xlabel('$w$'); ylabel('$|\mathbf{\hat{e}}_h|/|\mathbf{\hat{u}}|$')

width = 0;
h = 0;
ps = get(gcf, 'Position');
ratio = 1;
paperWidth = 10;
paperHeight = paperWidth*ratio;
set(gcf, 'paperunits', 'centimeters');
set(gcf, 'papersize', [paperWidth-width paperHeight-h]);
set(gcf, 'PaperPosition', [-width    -h   paperWidth+width paperHeight+h]);
print(gcf, '-dpdf', '../../tex/thesis/fig/wStability2.pdf');


%%
% r = .1;
% inside = sum(G.nodes.coords.^2,2) < r^2;
% 
% err = abs(sol1.nodeValues-uChi);
% err(inside) = max(err(~inside));
% sol = sol1;
% sol.nodeValues = err;
% sol = calculateCellAverages(G,sol);
% fig4 = figure;
% plotCellData(G, sol.cellMoments)
% colorbar;
% 
% err = abs(sol2.nodeValues-uChi);
% err(inside) = max(err(~inside));
% sol = sol2;
% sol.nodeValues = err;
% sol = calculateCellAverages(G,sol);
% fig5 = figure;
% plotCellData(G, sol.cellMoments)
% colorbar;

%%