function G = unitCubeTetrahedrons(dims)

nPx = dims(1)+1;
nPy = dims(2)+1;
nPz = dims(3)+1;

dx = 1/(nPx-1);
xx = 0:dx:1;
dy = 1/(nPy-1);
yy = 0:dy:1;
dz = 1/(nPz-1);
zz = 0:dz:1;

[x, y, z] = meshgrid(xx, yy, zz);   

x(:,2:nPx-1,:) = x(:,2:nPx-1,:) + random('Normal', 0, dx/4, nPy, nPx-2, nPz);
y(2:nPy-1,:,:) = y(2:nPy-1,:,:) + random('Normal', 0, dy/4, nPy-2, nPx, nPz);
z(:,:,2:nPz-1) = z(:,:,2:nPz-1) + random('Normal', 0, dz/4, nPy, nPx, nPz-2);


P = [x(:), y(:), z(:)];
T = DelaunayTri(P);

G = tetrahedralGrid(P, T.Triangulation);

% f = @(X) X(:,1);
% 
% G = computeGeometry(G);
% G = mrstGridWithFullMappings(G);
% G = computeVEMGeometry(G,f);
% 
% for i = 1:G.faces.num
%     nodeNum = G.faces.nodePos(i):G.faces.nodePos(i+1)-1;
%     nodes = G.faces.nodes(nodeNum);
%     X = G.nodes.coords(nodes,:)
%     plot3([X(:,1); X(1,1)], [X(:,2); X(1,2)],[X(:,3); X(1,3)], '-')
%     edgeNum = G.faces.edgePos(i):G.faces.edgePos(i+1)-1;
%     edges = G.faces.edges(edgeNum);
%     G.edges.lengths(edges)
% end

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