clc; close all;

run('../../matlab/project-mechanics-fractures/mystartup.m')

n = 10;
nx = n; ny = n;
G = unitSquare(nx, ny);
% G = cartGrid([nx,ny],[1,1]);
G = sortEdges(G);
G = mrstGridWithFullMappings(G);
G = computeGeometry(G);

% delta = 0.001;
% gD = @(X) 10*ones(size(X,1),1);
% gN = @(X) zeros(size(X,1),1);
% f = @(X) 1/sqrt(2*delta*pi).*exp(-((X(:,1)-0.5).^2 + (X(:,2)-0.5).^2)./(2*delta));
% 
% boundaryEdges = find((G.faces.neighbors(:,1) == 0) + (G.faces.neighbors(:,2) == 0));
% bNeu = zeros(numel(boundaryEdges),1);
% for e = 1:numel(boundaryEdges)
%     X = G.nodes.coords(G.faces.nodes(G.faces.nodePos(e):G.faces.nodePos(e +1)-1));
%     if X(:,1) < 0.001
%         bNeu(e) = 1;
%     end
% end
% bc = struct('bcFunc', {{gN,gD}}, 'bcFaces', {{boundaryEdges(bNeu == 1), boundaryEdges(bNeu == 0)}}, 'bcType', {{'neu', 'dir'}});

% gD = @(X) 1/4.*((X(:,1)-0.5).^2 + (X(:,2)-0.5).^2);
% f = @(X) -ones(size(X,1),1);
% boundaryEdges = find((G.faces.neighbors(:,1) == 0) + (G.faces.neighbors(:,2) == 0));
% bc = struct('bcFunc', {{gD}}, 'bcFaces', {{boundaryEdges}}, 'bcType', {{'dir'}});


% g = @(X) 1/6.*((X(:,1)-0.5).^3 + (X(:,2)-0.5).^3);
% f = @(X) -(X(:,1)-0.5) - (X(:,2)-0.5);

% g = @(X) X(:,2).*(1-X(:,2)).*X(:,1).^3;
% f = @(X) -6.*X(:,1).*X(:,2).*(1-X(:,2)) + 2.* X(:,1).^3;

gD = @(X) -log(1./((X(:,1)+0.1).^2 + (X(:,2)+0.1).^2));
f = @(X) zeros(size(X,1),1);
boundaryEdges = find((G.faces.neighbors(:,1) == 0) + (G.faces.neighbors(:,2) == 0));
bc = struct('bcFunc', {{gD}}, 'bcFaces', {{boundaryEdges}}, 'bcType', {{'dir'}});

U = VEM2D(G,f,bc);

h = 0;
Nc = G.cells.num;
baricenters = zeros(Nc,2);
for c = 1:Nc
    nodeNum = G.cells.nodePos(c):G.cells.nodePos(c+1)-1;        
    nodes = G.cells.nodes(nodeNum);
    X = G.nodes.coords(nodes,:);
    [~, baricenters(c,:)] = baric(X);
    h = max(h,cellDiameter(X));
end

X = [G.nodes.coords ; G.faces.centroids ; baricenters];
Uexact = gD(X);

errVec = U - Uexact;
err = sqrt(h)*norm(errVec, 2);
fprintf('Error: %e', err);

fig1 = figure;
plotGridWithDofs(G,bc);
fig2 = figure;
plotVEM(G, U, '');
fig3 = figure;
plotVEM(G, U, 'dof');

fig4 = figure;
plot(errVec);

%   Make function for computing baricenters of all cells
%   Make function for detecting boundary dofs
%   Change VEM2D to accept force term f