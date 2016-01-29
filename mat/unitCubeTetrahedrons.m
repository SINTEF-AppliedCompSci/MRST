function G = unitCubeTetrahedrons(n)

nP = n+1;

d = 1/(nP-1);
coords = 0:1/(nP-1):1;
[x, y, z] = meshgrid(coords, coords, coords);

% x(:,2:nP-1,:) = x(:,2:nP-1,:) + random('Normal', 0, d/4, nP, nP-2, nP);
% y(2:nP-1,:,:) = y(2:nP-1,:,:) + random('Normal', 0, d/4, nP-2, nP, nP);
% z(:,:,2:nP-1) = z(:,:,2:nP-1) + random('Normal', 0, d/4, nP, nP, nP-2);


P = [x(:), y(:), z(:)];
T = DelaunayTri(P);

G = tetrahedralGrid(P, T.Triangulation);

f = @(X) X(:,1);

G = computeGeometry(G);
G = mrstGridWithFullMappings(G);
G = computeVEMGeometry(G,f);

for i = 1:G.faces.num
    nodeNum = G.faces.nodePos(i):G.faces.nodePos(i+1)-1;
    nodes = G.faces.nodes(nodeNum);
    X = G.nodes.coords(nodes,:)
    plot3([X(:,1); X(1,1)], [X(:,2); X(1,2)],[X(:,3); X(1,3)], '-')
    edgeNum = G.faces.edgePos(i):G.faces.edgePos(i+1)-1;
    edges = G.faces.edges(edgeNum);
    G.edges.lengths(edges)
end

% plotGrid(G);
% hold on;
% plot3(G.edges.centroids(:,1), G.edges.centroids(:,2), G.edges.centroids(:,3),'o');
% plot3(G.faces.centroids(:,1), G.faces.centroids(:,2), G.faces.centroids(:,3),'o');
% add = .3;
% axis([-add, 1+add, -add, 1+add, -add, 1+add]); 

% for i = 1:G.cells.num;
%     faceNum = G.cells.facePos(i):G.cells.facePos(i+1)-1;
%     faces = G.cells.faces(faceNum);
%     nVec = bsxfun(@times, G.faces.normals(faces,:),...
%           ((-ones(size(faces,1),1)).^(G.faces.neighbors(faces,1) ~= i)));
%     
%     
%     plot3(G.faces.centroids(faces,1) + nVec(:,1), ...
%     G.faces.centroids(faces,2) + nVec(:,2), ...
%     G.faces.centroids(faces,3) + nVec(:,3), '*r');
% 
% end

% for i = 1:G.faces.num
% 
% nodeNum = G.faces.nodePos(i):G.faces.nodePos(i+1)-1;
% nodes = G.faces.nodes(nodeNum);
% X = G.nodes.coords(nodes,:);
% plot3([X(:,1);X(1,1)],[X(:,2);X(1,2)],[X(:,3);X(1,3)]);
% hold on;
% edgeNum = G.faces.edgePos(i):G.faces.edgePos(i+1)-1;
% edges = G.faces.edges(edgeNum);
% edgeN = G.faces.edgeNormals(edgeNum,:);
% nE = size(edges,1);
% for j = 1:nE
% plot3([G.edges.centroids(edges(j),1), G.edges.centroids(edges(j),1) + edgeN(j,1)], ...
%       [G.edges.centroids(edges(j),2), G.edges.centroids(edges(j),2) + edgeN(j,2)], ...
%       [G.edges.centroids(edges(j),3), G.edges.centroids(edges(j),3) + edgeN(j,3)], ...
%        '-r');
%    
% end
% axis equal
% hold off

end