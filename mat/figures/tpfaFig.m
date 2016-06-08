clc; clear all; close all;

addpath('../VEM2D/')
n = 5;
G = unitSquare([n,n],[10,10]);
G = computeGeometry(G);
G = mrstGridWithFullMappings(G);

dx = 10/(n-1);
dy = 10/(n-1);

addFac = .3;

intFaces = find(G.faces.centroids(:,1) > dx &  G.faces.centroids(:,1) < 10-dx & ...
                G.faces.centroids(:,2) > dy &  G.faces.centroids(:,2) < 10-dy);
% F = intFaces( max(G.faces.areas(intFaces)) == G.faces.areas(intFaces) );

nF = numel(intFaces);
F = intFaces(round(rand(1,1)*(nF-1) + 1));
Fc = G.faces.centroids(F,:);
Fn = G.faces.normals(F,:);
Fn = Fn/norm(Fn);
K = G.faces.neighbors(F,:)';
Kc = G.cells.centroids(K,:);

%%

figure;
plotGrid(G, K, 'faceAlpha', .2);

%%

hold on
plot([Kc(:,1);Fc(1)], [Kc(:,2);Fc(2)], 'ok', 'markerfacecolor', 'r', 'markersize', 6);
% quiver(Fc(1), Fc(2), Fn(1), Fn(2), 0,'k')
% quiver(Kc(1,1), Kc(1,2), Fc(1)-Kc(1,1), Fc(2)-Kc(1,2),0, 'k')

%%

fontsize = 8;

% set(gcf, 'DefaultTextInterpreter', 'LaTex');
% text(Kc(1,1), Kc(1,2),'$K_i$', 'fontsize', fontsize)
% text(Kc(2,1), Kc(2,2),'$K_j$', 'fontsize', fontsize)
% text(Kc(1,1), Kc(1,2),'$p_i$', 'fontsize', fontsize)
% text(Kc(2,1), Kc(2,2),'$p_j$', 'fontsize', fontsize)
% text((Kc(1,1)+Fc(1))/2, (Kc(1,2)+Fc(2))/2, '$\textbf{c}_{i,j}$', 'fontsize', fontsize)
% text((2*Fc(1)+Fn(1))/2, (2*Fc(2)+Fn(2))/2, '$\textbf{n}_{i,j}$', 'fontsize', fontsize)
% text(Fc(1), Fc(2), '$\pi_{i,j}$', 'fontsize', fontsize)
% text(Fc(1), Fc(2), '$F_{i,j}$', 'fontsize', fontsize)
axis equal off;

%%

% w = 0;
% h = 3;
% ps = get(gcf, 'Position');
% ratio = 1;
% paperWidth = 10;
% paperHeight = paperWidth*ratio;
% set(gcf, 'paperunits', 'centimeters');
% set(gcf, 'papersize', [paperWidth-w paperHeight-h]);
% set(gcf, 'PaperPosition', [-w    -h   paperWidth+w paperHeight+h]);
% print(gcf, '-dpdf', '../../tex/thesis/fig/tpfa_new.pdf');