clc; clear all; close all;
run('../../matlab/project-mechanics-fractures/mystartup.m')

nx = 30; ny = 30;
G = unitSquareV2(nx, ny);

G = mrstGridWithFullMappings(G);
G = computeGeometry(G);
Nc = G.cells.num;
Ne = G.faces.num;
Nn = G.nodes.num;
Ndof = Nn + Ne + Nc;

S = globS(G);

boundaryNodes = zeros(Ndof,1);
neighbors = G.faces.neighbors;
nodes = G.faces.nodes;
for e = 1:Ne
    if neighbors(e,1) == 0 || neighbors(e,2) == 0
        nodeNum = G.faces.nodePos(e):G.faces.nodePos(e+1)-1;
        boundaryNodes([G.faces.nodes(nodeNum)', Nn + e]) = 1;
    end
end

boundaryNodes = find(boundaryNodes);


SBC = zeros(numel(boundaryNodes),Ndof);
for i = 1:numel(boundaryNodes)
    SBC(i,boundaryNodes(i)) = 1;
end

SBC = eye(Ndof);

%S = S(activeNodes, activeNodes);

S(boundaryNodes,:) = SBC(boundaryNodes,:);

% g = @(X) -log(1./((X(:,1)+0.2).^2 + (X(:,2)+0.2).^2));
g = @(X) sin(X(:,1)*5*pi).*sin(X(:,2)*5*pi); %

baricenters = zeros(Nc,2);
for c = 1:Nc
    nodeNum = G.cells.nodePos(c):G.cells.nodePos(c+1)-1;        
    nodes = G.cells.nodes(nodeNum);
    X = G.nodes.coords(nodes,:);
    [~, baricenters(c,:)] = baric(X);
end


X = [G.nodes.coords ; G.faces.centroids ; baricenters];%G.cells.centroids];
b = zeros(Ndof, 1);  %ones(numel(activeNodes),1);
b(boundaryNodes) = g(X(boundaryNodes,:));

xi = S\b;

[x, y] = meshgrid(0:1/(2*nx-1):1, 0:1/(2*ny+1):1);
z = griddata(X(:,1), X(:,2), xi, x, y);

fig1 = figure;
plotGrid(G)
hold on
    N = X(boundaryNodes,:);
    plot(N(:,1), N(:,2),'or')
    plot(X(:,1),X(:,2),'.k')
hold off
fig2 = figure;
plotVEM(G, xi, 'dof');


xiExact = g(X);

norm(xi-xiExact)
