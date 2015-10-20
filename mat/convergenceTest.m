clc; close all;

Nk = 7;
errVec = zeros(Nk,1);
hVec = zeros(Nk,1);
for k = 1:Nk
n = 2^((k+6)/2);
G = cartGrid([n, n],[1,1]);
G = mrstGridWithFullMappings(G);
G = computeGeometry(G);

% g = @(X) zeros(size(X,1),1);
% f = @(X) exp(-((X(:,1)-0.5).^2 + (X(:,2)-0.5).^2)./0.01);
 
% g = @(X) 1/4.*((X(:,1)-0.5).^2 + (X(:,2)-0.5).^2);
% f = @(X) -ones(size(X,1),1);

g = @(X) X(:,2).*(1-X(:,2)).*X(:,1).^3;
f = @(X) -6.*X(:,1).*X(:,2).*(1-X(:,2)) + 2.* X(:,1).^3;

% g = @(X) -log(1./((X(:,1)+0.1).^2 + (X(:,2)+0.1).^2));
% f = @(X) zeros(size(X,1),1);

[U, hG] = VEM2D(G,f,g);

Nc = G.cells.num;
baricenters = zeros(Nc,2);
for c = 1:Nc
    nodeNum = G.cells.nodePos(c):G.cells.nodePos(c+1)-1;        
    nodes = G.cells.nodes(nodeNum);
    X = G.nodes.coords(nodes,:);
    [~, baricenters(c,:)] = baric(X);
end

X = [G.nodes.coords ; G.faces.centroids ; baricenters];
Uexact = g(X);

errVec(k) = norm(U - Uexact, Inf);
hVec(k) = hG;

end
loglog(hVec, errVec);

polyfit(log(hVec),log(errVec),1)
