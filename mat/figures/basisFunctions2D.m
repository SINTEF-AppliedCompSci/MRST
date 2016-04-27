clc; clear all; close all;

addpath('../VEM2D/')
addpath('../')

%%

G = unitSquare([5,5],[1,1]);

G = sortEdges(G);
G = computeVEM2DGeometry(G);
cut = 3;

bEdg = find(any(G.faces.neighbors == 0,2));


bc = VEM2D_addBC([], G, bEdg, 'pressure', 0);
[sol, G] = VEM2D(G, 0, 1, bc, 'projectors', true);

intCells = find(G.cells.centroids(:,1) > .2 & ...
                G.cells.centroids(:,1) < .8 & ...
                G.cells.centroids(:,2) > .2 & ...
                G.cells.centroids(:,2) < .8);
            
nK = numel(intCells);
K = intCells(round(rand(1,1)*(nK-1)+1));

PNstar = G.cells.PNstarT(G.cells.PNstarPos(K):G.cells.PNstarPos(K+1)-1,:)';
Kc = G.cells.centroids(K,:);
hK = G.cells.diameters(K);

[m, ~, ~] = retrieveMonomials(2,1);

nodes = G.cells.nodes(G.cells.nodePos(K):G.cells.nodePos(K+1)-1);
V = G.nodes.coords(nodes,:);
nN = numel(nodes);
xMax = max(V(:,1)); xMin = min(V(:,1));
yMax = max(V(:,2)); yMin = min(V(:,2));

[x,y] = meshgrid(linspace(xMin, xMax, 100), linspace(yMin, yMax, 100));
X = [x(:), y(:)];

X = X(inpolygon(X(:,1), X(:,2), V(:,1), V(:,2)),:);
X = bsxfun(@minus, X, Kc)/hK;
mVals = m(X);

for n = 1:nN
    PNphi = mVals*PNstar(:,n);
    plot3(X(:,1), X(:,2), PNphi)
end
%%

figure;

plotVEM2D(G, sol, 1, 'edgecolor', 'none');
view(-37.5, 40)
xlabel('$x$', 'interpreter', 'latex')
ylabel('$y$', 'interpreter', 'latex')
zlabel('$z$', 'interpreter', 'latex')
axis equal;
% axis([0 1 0 1 0 1])
% 
% set(gca,'XTick',[0 1] );
% set(gca,'YTick',[0 1] );
% set(gca,'ZTick',[0 1] );
    
%     ps = get(gcf, 'Position');
%     ratio = (ps(4)-ps(2)) / (ps(3)-ps(1));
%     paperWidth = 10;
%     paperHeight = paperWidth*ratio - cut;
%     set(gcf, 'paperunits', 'centimeters');
%     set(gcf, 'papersize', [paperWidth paperHeight]);
%     set(gcf, 'PaperPosition', [0    0   paperWidth paperHeight]);
%     dest = strcat('../../tex/thesis/fig/Phi', num2str(i), '.pdf');
%     print(gcf, '-dpdf', dest);
    


%%
figure;
fill(P(:,1), P(:,2), 'y', 'facealpha', .2)
hold on
xK = mean(polygonInt(G, 1:G.cells.num, @(X) X(:,1), 7)./G.cells.volumes);
yK = mean(polygonInt(G, 1:G.cells.num, @(X) X(:,2), 7)./G.cells.volumes);
plot(xK, yK, 'ok', 'MarkerFaceColor', 'r', 'markersize', 4);
text(xK*.9, yK*.9,'$\textbf{x}_K$', 'interpreter', 'latex')
plot([P(2,1),P(5,1)],[P(2,2),P(5,2)],'k--')
text(.65, .25,'$h_K$', 'interpreter', 'latex')
addFac = .1;
for i = 1:n
    cVec = P(i,:) - [xK, yK];
    x = P(i,1) + cVec(1)*addFac;
    y = P(i,2) + cVec(2)*addFac;
    text(x,y,num2str(i),'interpreter', 'latex', 'fontsize', 8);
end
axis equal off

%%

% ps = get(gcf, 'Position');
% ratio = (ps(4)-ps(2)) / (ps(3)-ps(1));
% paperWidth = 10;
% paperHeight = paperWidth*ratio - cut;
% set(gcf, 'paperunits', 'centimeters');
% set(gcf, 'papersize', [paperWidth paperHeight]);
% set(gcf, 'PaperPosition', [0    0   paperWidth paperHeight]);
% print(gcf, '-dpdf', '../../tex/thesis/fig/BasisElement.pdf');