clc; clear all; close all;

addpath('../VEM2D/')
addpath('../')

load('basisGrid3D');
%%


faces        = G.cells.faces(G.cells.facePos(K):G.cells.facePos(K+1)-1);
faceNormals  = G.faces.normals(faces,:);
nF = numel(faces);
faceSigns    = (-ones(nF,1)).^(G.faces.neighbors(faces,1) ~= K);
faceNormals  = bsxfun(@times, faceNormals,faceSigns);
faces        = faces(sum(bsxfun(@times, faceNormals, [0,0,1]), 2) < 0);

edges = G.cells.edges(G.cells.edgePos(K):G.cells.edgePos(K+1)-1);
allEdges = edges;
faceEdges    = G.faces.edges(mcolon(G.faces.edgePos(faces),G.faces.edgePos(faces+1)-1));
edges = edges(ismember(edges,faceEdges));
nodes = G.cells.nodes(G.cells.nodePos(K):G.cells.nodePos(K+1)-1);
faceNodes    = G.faces.nodes(mcolon(G.faces.nodePos(faces),G.faces.nodePos(faces+1)-1));
nodes = nodes(ismember(nodes,faceNodes));

Kc = G.cells.centroids(K,:);

X = G.nodes.coords(nodes,:);
Ec = G.edges.centroids(edges,:);
Fc = G.faces.centroids(faces,:);

nN = numel(nodes);
nE = numel(edges);
nF = numel(faces);

plotGrid(G, K, 'facealpha', .2, 'edgecolor', 'none')
hold on;

for i = 1:nE
    edgeNodes = G.edges.nodes(G.edges.nodePos(edges(i)):G.edges.nodePos(edges(i)+1)-1);
    Xe = G.nodes.coords(edgeNodes,:);
    plot3(Xe(:,1), Xe(:,2), Xe(:,3),'k');
end

X  = G.nodes.coords(G.edges.nodes(mcolon(G.edges.nodePos(edges),G.edges.nodePos(edges+1)-1)),:);
h1 = plot3(X(:,1), X(:,2), X(:,3),'ok','MarkerFaceColor', [0 0.4470 0.7410], 'markersize', 6);
Ec = G.edges.centroids(edges,:);
h2 = plot3(Ec(:,1),Ec(:,2), Ec(:,3),'sk', 'MarkerFaceColor', [0.8500 0.3250 0.0980], 'markersize', 6);
Fc = G.faces.centroids(faces,:);
h3 = plot3(Fc(:,1), Fc(:,2), Fc(:,3), 'dk', 'MarkerFaceColor', 'w', 'markersize', 6);
Kc = G.cells.centroids(K,:);
h4 = plot3(Kc(:,1), Kc(:,2), Kc(:,3), 'pk', 'MarkerFaceColor', 'm', 'markersize', 6);

remEdges = allEdges(~ismember(allEdges,edges));
nRE = numel(remEdges);
for i = 1:nRE
    edgeNodes = G.edges.nodes(G.edges.nodePos(remEdges(i)):G.edges.nodePos(remEdges(i)+1)-1);
    Xe = G.nodes.coords(edgeNodes,:);
    h = plot3(Xe(:,1), Xe(:,2), Xe(:,3));
    set(h, 'color', [0.5 0.5 0.5])
end

axis equal off
hold off



h = legend([h1, h2, h3, h4], '$\mathcal{V}^K$', '$\mathcal{E}^K$',...
                             '$\mathcal{F}^K$', '$\mathcal{P}^K$');
set(h, 'interpreter', 'latex')
set(h, 'fontsize', 14);
axis equal off;

%%

w = 0;
h = 2;
ps = get(gcf, 'Position');
ratio = 1;
paperWidth = 12;
paperHeight = paperWidth*ratio;
set(gcf, 'paperunits', 'centimeters');
set(gcf, 'papersize', [paperWidth-w paperHeight-h]);
set(gcf, 'PaperPosition', [-w    -h   paperWidth+w paperHeight+h]);
print(gcf, '-dpdf', '../../tex/thesis/fig/BasisElementDofs3D.pdf');