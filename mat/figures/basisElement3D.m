clc; clear all; close all;

addpath('../../../pebiGridding/voronoi3D/')
addpath('../VEM3D')
addpath('../')

G = voronoiCube(50,[1,1,1]);

G = computeVEM3DGeometry(G);

intCells = find(G.cells.centroids(:,1) > .2 & ...
                G.cells.centroids(:,1) < .8 & ...
                G.cells.centroids(:,2) > .2 & ...
                G.cells.centroids(:,2) < .8 & ...
                G.cells.centroids(:,3) > .2 & ...
                G.cells.centroids(:,3) < .8);
            
nK = numel(intCells);
K = intCells(round(rand(1,1)*(nK-1)+1));

% plotGrid(G, K, 'facealpha', .2)
% hold on
% Kc = G.cells.centroids(K,:);
% xK = Kc(1); yK = Kc(2); zK = Kc(3);
% plot3(xK, yK, zK, 'ok', 'MarkerFaceColor', 'r', 'markersize', 4);
% text(xK*.98, yK*.98, zK*.98, '$\textbf{x}_K$', 'interpreter', 'latex')
% axis equal off
% hold off


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

remEdges = allEdges(~ismember(allEdges,edges));
nRE = numel(remEdges);
for i = 1:nRE
    edgeNodes = G.edges.nodes(G.edges.nodePos(remEdges(i)):G.edges.nodePos(remEdges(i)+1)-1);
    Xe = G.nodes.coords(edgeNodes,:);
    h = plot3(Xe(:,1), Xe(:,2), Xe(:,3));
    set(h, 'color', [0.5 0.5 0.5])
end

xK = Kc(1); yK = Kc(2); zK = Kc(3);
plot3(xK, yK, zK, 'ok', 'MarkerFaceColor', 'r', 'markersize', 4);
text(xK*.98, yK*.98, zK*.98, '$\textbf{x}_K$', 'interpreter', 'latex')
axis equal off
hold off

%%

cut = 4;
ps = get(gcf, 'Position');
ratio = (ps(4)-ps(2)) / (ps(3)-ps(1));
paperWidth = 10;
paperHeight = paperWidth*ratio - cut;
set(gcf, 'paperunits', 'centimeters');
set(gcf, 'papersize', [paperWidth paperHeight]);
set(gcf, 'PaperPosition', [0    0   paperWidth paperHeight]);
print(gcf, '-dpdf', '../../tex/thesis/fig/BasisElement3D_new.pdf', '-r0');

%%

G = computeVEM3DGeometry(G);

bF = find(any(G.faces.neighbors == 0,2));
bc = VEM3D_addBC([], bF, 'pressure', 0);

[sol, G] = VEM3D(G, 0, bc, 2, 'faceProjectors', true);

%%


faces = G.cells.faces(G.cells.facePos(K):G.cells.facePos(K+1)-1);
faceNormals = G.faces.normals(faces,:);

faces = faces(dot(G.faces.normals(faces,:),[0,0,1])) > 0;
nF = numel(faces);

nodes = G.cells.nodes(G.cells.nodePos(K):G.cells.nodePos(K+1)-1);
nodes = nodes(G.nodes.coords(nodes,3) < G.cells.centroids(K,3)*.8);
nN = numel(nodes);
figure;
m = retrieveMonomials(2,2);
for k = 1:nN
    clf;
    n = nodes(k);
    plotGrid(G,K, 'facecolor', 'none')
    for i = 1:nF
            F = faces(i);
            PNstar = G.faces.PNstarT(G.faces.PNstarPos(F):G.faces.PNstarPos(F+1)-1,:)';
            T = G.faces.localCoords(G.faces.TPos(F):G.faces.TPos(F+1)-1,:);
            b = G.faces.centroids(F,:);
            hF = G.faces.diameters(F);
            faceNodes = G.faces.nodes(G.faces.nodePos(F):G.faces.nodePos(F+1)-1);
            V = G.nodes.coords(faceNodes,:);
            V = ((bsxfun(@minus, V,b))*T);
            np = hF*100;
            x = linspace(min(V(:,1)), max(V(:,1)), np);
            y = linspace(min(V(:,2)), max(V(:,2)), np);
            [x,y] = meshgrid(x,y);
            X = [x(:), y(:)];
            X = [X(inpolygon(X(:,1), X(:,2), V(:,1), V(:,2)),:)];
            nNF = numel(faceNodes);
            for j = 1:nNF
                P1 = V(j,:);
                P2 = V(mod(j,nNF)+1,:); 
                npn = norm(P1-P2,2)*100;
                alpha = linspace(0,1,npn)';
                VV = [P2(1)*alpha,     P2(2)*alpha] ...
                   + [P1(1)*(1-alpha), P1(2)*(1-alpha)      ];
               X = [X;VV];
            end
            tri = delaunay(X);
            Xmon = X/hF;
            nn = find(faceNodes == n);
            if ~isempty(nn)
                PNphi = m(Xmon)*PNstar(:,nn);
            else
                PNphi = m(Xmon)*zeros(size(PNstar,1),1);
            end
            
            X = bsxfun(@plus, X*T', b);
            hold on
            trisurf(tri,X(:,1), X(:,2), X(:,3), PNphi, 'edgecolor', 'none')
    end
    axis equal off;
    colorbar;
    cut = 4;
    ps = get(gcf, 'Position');
    ratio = (ps(4)-ps(2)) / (ps(3)-ps(1));
    paperWidth = 10;
    paperHeight = paperWidth*ratio - cut;
    set(gcf, 'paperunits', 'centimeters');
    set(gcf, 'papersize', [paperWidth paperHeight]);
    set(gcf, 'PaperPosition', [0    0   paperWidth paperHeight]);
    dest = strcat('../../tex/thesis/fig/basis3D/Phi', num2str(k),'3D_new.pdf');
    print(gcf, '-dpdf', dest, '-r0');
end