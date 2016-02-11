function p = plotCellDetails(G,K)

close all

Kc = G.cells.centroids(K,:);

nodeNum = G.cells.nodePos(K):G.cells.nodePos(K+1)-1;
nodes   = G.cells.nodes(nodeNum);
X       = G.nodes.coords(nodes,:);

edgeNum = G.cells.edgePos(K):G.cells.edgePos(K+1)-1;
edges   = G.cells.edges(edgeNum,:);
Ec      = G.edges.centroids(edges,:);

faceNum = G.cells.facePos(K):G.cells.facePos(K+1)-1;
faces   = G.cells.faces(faceNum);
nF = size(faces,1);
Fc      = G.faces.centroids(faces,:);
faceNormals = G.faces.normals(faces,:);
faceSigns    = (-ones(nF,1)).^(G.faces.neighbors(faces,1) ~= K);
faceNormals = bsxfun(@times, faceNormals,faceSigns);
edgeNum = mcolon(G.faces.edgePos(faces),G.faces.edgePos(faces+1)-1);
edgeNormals = G.faces.edgeNormals(edgeNum,:);


plotGrid(G,K,'FaceAlpha', .5);
hold on
plot3(X(:,1), X(:,2), X(:,3), '.b', 'MarkerSize', 12);
plot3(Kc(1), Kc(2), Kc(3), '.r', 'MarkerSize', 13);
for i = 1:nF
    plot3([Fc(i,1), Fc(i,1) + faceNormals(i,1)], ...
          [Fc(i,2), Fc(i,2) + faceNormals(i,2)], ...
          [Fc(i,3), Fc(i,3) + faceNormals(i,3)], '-ok', 'MarkerSize', 5)
end
hold off

axis equal off

p = 1;

end