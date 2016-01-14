function [faces, midpoints, normals] = faceData(G,K)

faceNum = G.cells.facePos(K):G.cells.facePos(K+1)-1;
faces = G.cells.faces(faceNum);
if size(faces,1) == 1
    faces = faces';
end
midpoints = G.faces.centroids(faces,:);
                            %   Area wheighted edge normals.
normals = G.faces.normals(faces,:);
neighbors = G.faces.neighbors(faces,:);
                            %   Fix orientation of edge normals.
m = (-ones(length(normals),1)).^(neighbors(:,1) ~= K);
normals = [m,m,m].*normals;

end