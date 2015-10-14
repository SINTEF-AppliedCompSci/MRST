clc; clear all; close all;

run('../../matlab/project-mechanics-fractures/mystartup.m')

nx = 30; ny = 30;
G = unitSquare(nx, ny);
G = mrstGridWithFullMappings(G);
G = computeGeometry(G);

% g = @(X) zeros(size(X,1),1);
% f = @(X) exp(-((X(:,1)-0.5).^2 + (X(:,2)-0.5).^2));

g = @(X) 1/4.*(X(:,1).^2 + X(:,2).^2);
f = @(X) -ones(size(X,1),1);

% g = @(X) -log(1./((X(:,1)+0.1).^2 + (X(:,2)+0.1).^2));
% f = @(X) zeros(size(X,1),1);

U = VEM2D(G,f,g);

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

norm(U - Uexact)

fig1 = figure;
plotGridWithDofs(G);
fig2 = figure;
plotVEM(G, U, '');
fig3 = figure;
plotVEM(G, U, 'dof');
fig4 = figure;
plotVEM(G,Uexact, '');

%   Make function for computing baricenters of all cells
%   Make function for detecting boundary dofs
%   Change VEM2D to accept force term f