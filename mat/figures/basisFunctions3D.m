clc; clear all; close all;

addpath('../../../pebiGridding/voronoi3D/')
addpath('../VEM3D')
addpath('../VEM2D')
addpath('../')

% G = voronoiCube(50,[1,1,1]);
% 
% G = computeVEM3DGeometry(G);
% 
% intCells = find(G.cells.centroids(:,1) > .2 & ...
%                 G.cells.centroids(:,1) < .8 & ...
%                 G.cells.centroids(:,2) > .2 & ...
%                 G.cells.centroids(:,2) < .8 & ...
%                 G.cells.centroids(:,3) > .2 & ...
%                 G.cells.centroids(:,3) < .8);
%             
% nK = numel(intCells);
% K = intCells(round(rand(1,1)*(nK-1)+1));
% 
% % plotGrid(G, K, 'facealpha', .2)
% % hold on
% % Kc = G.cells.centroids(K,:);
% % xK = Kc(1); yK = Kc(2); zK = Kc(3);
% % plot3(xK, yK, zK, 'ok', 'MarkerFaceColor', 'r', 'markersize', 4);
% % text(xK*.98, yK*.98, zK*.98, '$\textbf{x}_K$', 'interpreter', 'latex')
% % axis equal off
% hold off

%%

load('basisGrid3D.mat')

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
text(xK*.98, yK*.95, zK*.98, '$\textbf{x}_K$', 'interpreter', 'latex')
axis equal off
hold off

%%
% 
% save('basisGrid3D.mat', 'G', 'K')

%%

% cut = 4;
% ps = get(gcf, 'Position');
% ratio = (ps(4)-ps(2)) / (ps(3)-ps(1));
% paperWidth = 10;
% paperHeight = paperWidth*ratio - cut;
% set(gcf, 'paperunits', 'centimeters');
% set(gcf, 'papersize', [paperWidth paperHeight]);
% set(gcf, 'PaperPosition', [0    0   paperWidth paperHeight]);
% print(gcf, '-dpdf', '../../tex/thesis/fig/basis3D/BasisElement3D_new.pdf');

%%

G = computeVEM3DGeometry(G);

bF = find(any(G.faces.neighbors == 0,2));
bc = VEM3D_addBC([], bF, 'pressure', 0);

[sol, G] = VEM3D(G, 0, bc, 2, 'faceProjectors', true);

%%

faces        = G.cells.faces(G.cells.facePos(K):G.cells.facePos(K+1)-1);
faceNormals  = G.faces.normals(faces,:);
nF = numel(faces);
faceSigns    = (-ones(nF,1)).^(G.faces.neighbors(faces,1) ~= K);
faceNormals  = bsxfun(@times, faceNormals,faceSigns);
faces        = faces(sum(bsxfun(@times, faceNormals, [0,0,1]), 2) < 0);
faceNodes    = G.faces.nodes(mcolon(G.faces.nodePos(faces),G.faces.nodePos(faces+1)-1));
nF = numel(faces);
nodes = G.cells.nodes(G.cells.nodePos(K):G.cells.nodePos(K+1)-1);
nodes = nodes(ismember(nodes,faceNodes));
nN = numel(nodes);
figure;
m = retrieveMonomials(2,2);
for k = 1:nN
    clf;
    n = nodes(k);
%     plotGrid(G,K, 'facecolor', 'none')
    for i = 1:nE
        edgeNodes = G.edges.nodes(G.edges.nodePos(edges(i)):G.edges.nodePos(edges(i)+1)-1);
        Xe = G.nodes.coords(edgeNodes,:);
        plot3(Xe(:,1), Xe(:,2), Xe(:,3),'k');
        hold on
    end

    for i = 1:nF
        
        F = faces(i);
        T = G.faces.localCoords(G.faces.TPos(F):G.faces.TPos(F+1)-1,:);
        b = G.faces.centroids(F,:);
        hF = G.faces.diameters(F);
        faceNodes = G.faces.nodes(G.faces.nodePos(F):G.faces.nodePos(F+1)-1);
        V = G.nodes.coords(faceNodes,:);

        nn = find(faceNodes == n);
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
        
%         for j = 1:nNF
%             xx = edgePts{j};
%             for l = 1:size(xx,1)
%                 remNodes = sum(bsxfun(@minus, X,xx(l,:)).^2,2) < 1/(4*np^2);
%                 X = X(~remNodes,:);
%                 X = [X;xx];
%             end
%         end
        tri = delaunay(X);
        H = triangleGrid(X,tri);
        H = computeGeometry(H);
        Kc = mean(H.cells.volumes);
        remCells = find(H.cells.volumes < .05*Kc);
        H = removeCells(H,remCells);
        H = computeVEM2DGeometry(H);
        
        gD = 0;
        if ~isempty(nn)

            vec = [V(nn,:)-V(mod(nn-2,nNF)+1,:); V(nn,:) - V(mod(nn,nNF)+1,:)];
            A = vec\[1;1];
            B = -V(mod(nn,nNF)+1,:)*A;
            gD = @(X) X*A + B;
        end
            
        bF = find(any(H.faces.neighbors == 0,2));
        
        bFdof = true(numel(bF),1);
        if ~isempty(nn)
            d1 = sqrt(sum(bsxfun(@minus, H.faces.centroids(bF,:),V(nn,:)).^2,2));
            d2 = sqrt(sum(bsxfun(@minus, H.faces.centroids(bF,:),V(mod(nn-2,nNF)+1,:)).^2,2));
            d3 = sqrt(sum(bsxfun(@minus, H.faces.centroids(bF,:),V(mod(nn,nNF)+1,:)).^2,2));
            dV1 = norm(vec(1,:),2);
            dV2 = norm(vec(2,:),2);

            tol = 1e-3;
            bFdof = any([abs(bsxfun(@minus, d1 + d2, dV1))./dV1, abs(bsxfun(@minus, d1 + d3, dV2))./dV2] < tol,2);
        end

        bc = VEM2D_addBC([], H, bF(bFdof), 'pressure', gD);
        bc = VEM2D_addBC(bc, H, bF(~bFdof), 'pressure', 0);

        sol = VEM2D(H,0,bc,2);

        X = H.nodes.coords;
        tri = delaunay(X);
        X = bsxfun(@plus, X*T', b);
        hold on
        trisurf(tri,X(:,1), X(:,2), X(:,3), sol.nodeValues, 'edgecolor', 'none')
%         figure;
%         plotGrid(H)
%         hold on
%         plot(H.faces.centroids(bF(bFdof),1), H.faces.centroids(bF(bFdof),2),'o');
%         plot(V(nn,1), V(nn,2),'o', 'markerfacecolor', 'r')
%         plot(H.faces.centroids(bF,1), H.faces.centroids(bF,2),'r.')
    end

%     axis equal off;
%     colorbar;
%     cut = 4;
%     ps = get(gcf, 'Position');
%     ratio = (ps(4)-ps(2)) / (ps(3)-ps(1));
%     paperWidth = 10;
%     paperHeight = paperWidth*ratio - cut;
%     set(gcf, 'paperunits', 'centimeters');
%     set(gcf, 'papersize', [paperWidth paperHeight]);
%     set(gcf, 'PaperPosition', [0    0   paperWidth paperHeight]);
%     dest = strcat('../../tex/thesis/fig/basis3D/Phi', num2str(k),'3D_new.pdf');
%     print(gcf, '-dpdf', dest, '-r0');
      for i = 1:nE
        edgeNodes = G.edges.nodes(G.edges.nodePos(edges(i)):G.edges.nodePos(edges(i)+1)-1);
        Xe = G.nodes.coords(edgeNodes,:);
        plot3(Xe(:,1), Xe(:,2), Xe(:,3),'k');
        hold on
      end

axis equal off;
view([0,90])
colorbar;

h = 1;
w = 0;
ps = get(gcf, 'Position');
ratio = 1;
paperWidth = 10;
paperHeight = paperWidth*ratio;
set(gcf, 'paperunits', 'centimeters');
set(gcf, 'papersize', [paperWidth-w paperHeight-h]);
set(gcf, 'PaperPosition', [-w    -h   paperWidth+w paperHeight+h]);
dest = strcat('../../tex/thesis/fig/basis3D/Phi3D_', num2str(k), '.pdf');

print(gcf, '-dpdf', dest, '-r10');

end