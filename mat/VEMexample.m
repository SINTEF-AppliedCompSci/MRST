%--------------------------------------------------------------------------
%   Example of VEM for the Poisson problem.
%--------------------------------------------------------------------------

clc; close all;

run('../../matlab/project-mechanics-fractures/mystartup.m')

                            %   Generate rando polygon grid of roughly
                            %   (n-1)x(n-1) cells.
n = 15;
nx = n; ny = n;
G = unitSquare(nx, ny);
G = sortEdges(G);
G = mrstGridWithFullMappings(G);
G = computeGeometry(G);

                            %   Define source term and Dirichelt and
                            %   Nueumann BC's.
f =  @(X)  4*pi^2*sin(X(:,1)*pi).*cos(X(:,2)*pi);
gD = @(X)  2*sin(X(:,1)*pi).*cos(X(:,2)*pi) - log(1./((X(:,1)+0.1).^2 + (X(:,2)+0.1).^2));
gN = @(X) -2*pi*cos(X(:,1)*pi).*cos(X(:,2)*pi) - 2*(X(:,1)+0.1)./((X(:,1)+0.1).^2 + (X(:,2)+0.1).^2);

                            %   Idefntify boundary edges.
tol = 1e-6;
boundaryEdges = find((G.faces.neighbors(:,1) == 0) + ...
                     (G.faces.neighbors(:,2) == 0));
bNeu = abs(G.faces.centroids(boundaryEdges,1)) < tol;
bc = struct('bcFunc', {{gN, gD}}, 'bcFaces', ...
            {{boundaryEdges(bNeu), boundaryEdges(~bNeu)}}, ...
            'bcType', {{'neu', 'dir'}});

                            %   Calculate VEM solution.
U = VEM2D(G,f,bc);

                            %   Plot grid.
figure();
Nc = G.cells.num;
plotGrid(G, 'facecolor', 'w');
xlabel('x');
ylabel('y');
axis equal off;

                            %   Plot solution.
figure();
plotVEM(G, U, '');
view(-45,38);
axis([0, 1, 0, 1, -4, 1.5]);
set(gcf,'Units','normalized')