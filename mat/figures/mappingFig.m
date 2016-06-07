clc; clear all; close all;

addpath('../VEM2D/')
n = 4;
G = unitSquare([n,n],[1,1]);
G = computeGeometry(G);
G = mrstGridWithFullMappings(G);

dx = 1/(n-1);
dy = 1/(n-1);

addFac = .3;

intCells = find(G.cells.centroids(:,1) > dx   & ...
                G.cells.centroids(:,1) < 1-dx & ...
                G.cells.centroids(:,2) > dy   & ...
                G.cells.centroids(:,2) < 1-dy       );                 

nK = numel(intCells);
K = intCells(round(rand(1,1)*(nK-1) + 1));
Kc = G.cells.centroids(K,:);

nodeNum = G.cells.nodePos(K):G.cells.nodePos(K+1)-1;
nodes   = G.cells.nodes(nodeNum);
X       = G.nodes.coords(nodes,:);
nN = numel(nodes);

%%

figure;
subplot(1,2,1);
plotGrid(G, K, 'facecolor', 'k', 'faceAlpha', .2);
hold on;
plot(X(:,1), X(:,2), 'ok', 'MarkerFaceColor', 'r', 'markersize', 4);
% text(Kc(1), Kc(2),'$K$', 'interpreter', 'latex')
% for i = 1:nN
%     cVec = X(i,:) - Kc;
%     x = X(i,1) + cVec(1)*addFac*.9;
%     y = X(i,2) + cVec(2)*addFac*.9;
%     text(x,y,num2str(i),'interpreter', 'latex', 'fontsize', 5);
% end
axis equal off;

subplot(1,2,2);
plotGrid(G, 'faceAlpha', .2);
hold on;
plotGrid(G, K, 'facecolor', 'w');
plotGrid(G, K, 'facecolor', 'k', 'faceAlpha', .2);
plot(X(:,1), X(:,2), 'ok', 'MarkerFaceColor', 'r', 'markersize', 4);
% text(Kc(1), Kc(2),'$K$', 'interpreter', 'latex')
% for i = 1:nN
%     cVec = X(i,:) - Kc;
%     x = X(i,1) + cVec(1)*addFac;
%     y = X(i,2) + cVec(2)*addFac;
%     text(x,y,num2str(nodes(i)),'interpreter', 'latex', 'fontsize', 5);
% end
axis equal off;

% annotation(gcf,'arrow',[.46 .56], ...
%                            [.5  .5]);
% text(-.21, .5, '$\delta$', 'interpreter', 'latex')


cut = 4;

%%

w = 0;
h = 5;
ps = get(gcf, 'Position');
ratio = 1;
paperWidth = 10;
paperHeight = paperWidth*ratio;
set(gcf, 'paperunits', 'centimeters');
set(gcf, 'papersize', [paperWidth-w paperHeight-h]);
set(gcf, 'PaperPosition', [-w    -h   paperWidth+w paperHeight+h]);
print(gcf, '-dpdf', '../../tex/thesis/fig/Mapping_new.pdf');

%%

% save('basisElement.mat', 'G', 'K') 