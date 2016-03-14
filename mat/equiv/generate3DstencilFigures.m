clc; clear; close all;

G = cartGrid([2,2,2]);
k = 1;
f = @(X) X(:,1);
add = 1;
lineWidth = 2;

X = G.nodes.coords;
X = bsxfun(@minus, X,[1,1,1]);
G.nodes.coords = X;

G = computeVEMGeometry(G,f,k);

edgeNodes = G.edges.nodes;
nN = size(edgeNodes,1)/2;
Xe = G.nodes.coords(edgeNodes,:);
Xe = reshape(Xe,2,[]);


fig1 = figure();
plotGrid(G,'FaceAlpha', .3);
edgeNodes = G.edges.nodes;
nN = size(edgeNodes,1)/2;
hold on
plot3(Xe(:,1:nN), Xe(:,nN+1:2*nN), Xe(:,2*nN+1:3*nN),'k')
X1 = X((X(:,1)).^2 + (X(:,2)).^2 + (X(:,3)).^2 <= 1,:);
plot3(X1(:,1), X1(:,2), X1(:,3), 'ok', 'MarkerFaceColor', 'r')
line([-(1+add) 1+add], [0,0], [0,0], 'LineWidth', lineWidth);
line([0,0], [-(1+add) 1+add], [0,0], 'LineWidth', lineWidth);
line([0,0], [0,0], [-(1+add) 1+add], 'LineWidth', lineWidth);
axis equal off
view(3)

fig2 = figure();
plotGrid(G,'FaceAlpha', .3);
edgeNodes = G.edges.nodes;
nN = size(edgeNodes,1)/2;
Xe = G.nodes.coords(edgeNodes,:);
Xe = reshape(Xe,2,[]);
hold on
plot3(Xe(:,1:nN), Xe(:,nN+1:2*nN), Xe(:,2*nN+1:3*nN),'k')
X2 = X( ((X(:,1).^2 + X(:,3).^2 <= 2) & X(:,2) == 0) | ...
        ((X(:,2).^2 + X(:,3).^2 <= 2) & (X(:,1) == 0)| ...
        ((X(:,1).^2 + X(:,2).^2 <= 2) & (X(:,3) == 0)) ),:);
plot3(X2(:,1), X2(:,2), X2(:,3), 'ok', 'MarkerFaceColor', 'r')
line([-(1+add) 1+add], [-(1+add) 1+add], [0,0], 'LineWidth', lineWidth);
line([-(1+add) 1+add], [(1+add) -(1+add)], [0,0], 'LineWidth', lineWidth);
line([-(1+add) 1+add], [0,0], [-(1+add) 1+add], 'LineWidth', lineWidth);
line([-(1+add) 1+add], [0,0], [(1+add) -(1+add)], 'LineWidth', lineWidth);
line([0,0],[-(1+add) 1+add], [-(1+add) 1+add], 'LineWidth', lineWidth);
line([0,0],[-(1+add) 1+add], [(1+add) -(1+add)], 'LineWidth', lineWidth);

axis equal off
view(3)

fig3 = figure();
plotGrid(G,'FaceAlpha', .3);
edgeNodes = G.edges.nodes;
nN = size(edgeNodes,1)/2;
Xe = G.nodes.coords(edgeNodes,:);
Xe = reshape(Xe,2,[]);
hold on
plot3(Xe(:,1:nN), Xe(:,nN+1:2*nN), Xe(:,2*nN+1:3*nN),'k')
X3 = [X( X(:,1).^2 + X(:,2).^2 + X(:,3).^2 == 3, :); [0,0,0]];
plot3(X3(:,1), X3(:,2), X3(:,3), 'ok', 'MarkerFaceColor', 'r')
line([-(1+add) 1+add], [-(1+add) 1+add], [-(1+add) 1+add], 'LineWidth', lineWidth);
line([(1+add) -(1+add)], [-(1+add) 1+add], [-(1+add) 1+add], 'LineWidth', lineWidth);
line([-(1+add) 1+add], [(1+add) -(1+add)], [-(1+add) 1+add], 'LineWidth', lineWidth);
line([-(1+add) (1+add)], [-(1+add) (1+add)], [(1+add) -(1+add)], 'LineWidth', lineWidth);
axis equal off
view(3)