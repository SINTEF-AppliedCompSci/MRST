clc; clear all; close all;
run('../../matlab/project-mechanics-fractures/mystartup.m')

nx =20; ny = 20;              %   Grid dimensions.
G = cartGrid([nx, ny]);

G = mrstGridWithFullMappings(G);
G = computeGeometry(G);
Nc = G.cells.num;
Ne = G.faces.num;
Nn = G.nodes.num;
Ndof = Nn + Ne + Nc;

S = globS(G, nx, ny);

boundaryNodes = ones(Ndof,1);
neighbors = G.faces.neighbors;
nodes = G.faces.nodes;
for e = 1:Ne
    if neighbors(e,1) == 0 || neighbors(e,2) == 0
        boundaryNodes([nodes([2*e-1,2*e])', Nn + e]) = 1;
    end
end
boundaryNodes = find(boundaryNodes);
A = zeros(numel(boundaryNodes),Ndof);
for i = 1:numel(boundaryNodes)
    A(i,boundaryNodes(i)) = 1;
end

%S = S(activeNodes, activeNodes);

S(boundaryNodes,:) = A;

g = @(X) -log(1./((X(:,1)+0.2).^2 + (X(:,2)+0.2).^2));
nodes = [G.nodes.coords ; G.faces.centroids ; G.cells.centroids];
b = zeros(numel(boundaryNodes), 1);  %ones(numel(activeNodes),1);
b(boundaryNodes) = g(nodes(boundaryNodes,:));

xi = g(nodes);

xi(boundaryNodes) = S\b;

X = [G.nodes.coords; G.faces.centroids; G.cells.centroids];
[x,y] = meshgrid(0:0.5:nx, 0:0.5:ny);
z = griddata(X(:,1), X(:,2), xi, x, y);
surf(x,y,z);

xiExact = g(X);

norm(xi-xiExact)
