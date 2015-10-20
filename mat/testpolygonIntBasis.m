clc; close all;

G = cartGrid([1,1]);
G = mrstGridWithFullMappings(G);
G = computeGeometry(G);


K = 1;
vol = G.cells.volumes(K);
nodeNum = G.cells.nodePos(K):G.cells.nodePos(K+1)-1;
nodes = G.cells.nodes(nodeNum);
edgeNum = G.cells.facePos(K):G.cells.facePos(K+1)-1;
edges = G.cells.faces(edgeNum);

X = G.nodes.coords(nodes,:);
Xmid = G.faces.centroids(edges,:);
[~, XB] = baric(X);
xK = XB(1); yK = XB(2);
Xdof = [X; Xmid; XB];
                            %   Element diameter.
hK = 0;
n = size(X,1);
for i = 1:n
    hK = max(norm(repmat(X(i,:),n,1)-X),hK);
end
m =      @(X) 1/hK.*[X(:,1)-xK, ...                             %   (1,0)
                     X(:,2)-yK];                             %   (0,1)

I = polygonIntBasis(X,Xmid,XB,vol,m)