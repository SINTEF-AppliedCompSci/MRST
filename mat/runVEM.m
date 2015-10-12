clc; clear all; close all;
run('../../matlab/project-mechanics-fractures/mystartup.m')

% nx =20; ny = 20;              %   Grid dimensions.
% G = cartGrid([nx, ny],[1,1]);
n = 10;
G = unitSquare(n);
nx = n; ny = n;

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
        boundaryNodes([nodes([2*e-1,2*e])', Nn + e]) = 1;
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

g = @(X) -log(1./((X(:,1)+0.2).^2 + (X(:,2)+0.2).^2));% sin(X(:,1)*5*pi) + sin(X(:,2)*5*pi); %
X = [G.nodes.coords ; G.faces.centroids ; G.cells.centroids];
b = zeros(Ndof, 1);  %ones(numel(activeNodes),1);
b(boundaryNodes) = g(X(boundaryNodes,:));

xi = S\b;

[x, y] = meshgrid(0:1/(2*nx-1):1, 0:1/(2*ny+1):1);
z = griddata(X(:,1), X(:,2), xi, x, y);

fig1 = figure;
plotGrid(G)
fig2 = figure;
surf(x,y,z)
fig3 = figure;
plotNodeData(G,xi);


xiExact = g(X);

norm(xi-xiExact)
