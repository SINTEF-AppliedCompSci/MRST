function [edges, midpoints, normals] = faceData(G,K)

faceNum = G.cells.facePos(K):G.cells.facePos(K+1)-1;
edges = G.cells.faces(faceNum);
midpoints = G.faces.centroids(edges,:);
                            %   Area wheighted edge normals.
normals = G.faces.normals(edges,:);
neighbors = G.faces.neighbors(edges,:);
                            %   Fix orientation of edge normals.
m = (-ones(length(normals),1)).^(neighbors(:,1) ~= K);
normals = [m,m].*normals;

end