clc; close all;

Nk = 4;
errVec = zeros(Nk,1);
hVec = zeros(Nk,1);
for k = 1:Nk
n = 2^(k+4);
% G = cartGrid([n, n],[1,1]);
G = unitSquare(floor(n),floor(n));
G = mrstGridWithFullMappings(G);
G = computeGeometry(G);

% gD = @(X) 1/4.*((X(:,1)-0.5).^2 + (X(:,2)-0.5).^2);
% f = @(X) -ones(size(X,1),1);
% boundaryEdges = find((G.faces.neighbors(:,1) == 0) + (G.faces.neighbors(:,2) == 0));
% bc = struct('bcFunc', {{gD}}, 'bcFaces', {{boundaryEdges}}, 'bcType', {{'dir'}});

gD = @(X) -log(1./((X(:,1)+1).^2 + (X(:,2)+1).^2));
f = @(X) zeros(size(X,1),1);
boundaryEdges = find((G.faces.neighbors(:,1) == 0) + (G.faces.neighbors(:,2) == 0));
bc = struct('bcFunc', {{gD}}, 'bcFaces', {{boundaryEdges}}, 'bcType', {{'dir'}});

% gD = @(X) -log(1./((X(:,1)+0.1).^2 + (X(:,2)+0.1).^2));
% f = @(X) zeros(size(X,1),1);
% boundaryEdges = find((G.faces.neighbors(:,1) == 0) + (G.faces.neighbors(:,2) == 0));
% bc = struct('bcFunc', {{gD}}, 'bcFaces', {{boundaryEdges}}, 'bcType', {{'dir'}});

[U, hG] = VEM2D(G,f,bc);

Nc = G.cells.num;
baricenters = zeros(Nc,2);
hM = 0;
for c = 1:Nc
    nodeNum = G.cells.nodePos(c):G.cells.nodePos(c+1)-1;        
    nodes = G.cells.nodes(nodeNum);
    X = G.nodes.coords(nodes,:);
    [~, baricenters(c,:)] = baric(X);
    hM = hM + cellDiameter(X);
end

hM = hM/Nc;

X = [G.nodes.coords ; G.faces.centroids ; baricenters];
Uexact = gD(X);

errVec(k) = hM*norm(U - Uexact, Inf);
hVec(k) = hM;

end
loglog(hVec, errVec);

polyfit(log(hVec),log(errVec),1)
