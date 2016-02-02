clc; clear; close all;

run('../../matlab/project-mechanics-fractures/mystartup.m')

n = 2;


G = cartGrid([n,n,n],[2,2,2]);
% G = unitCubeTetrahedrons(n);

% %--------------------------------------------------------------------------
% %   -\delta u = 0,
% %           u = 1/(2\pi||x-C||)
% %--------------------------------------------------------------------------
% f = @(X) zeros(size(X,1),1);
% C = -[.2,.2,.2];
% gD = @(X) -1./(2*pi*sqrt(sum((X-repmat(C,size(X,1),1)).^2,2)));

%--------------------------------------------------------------------------
%   -\delta u = 1,
%           u = -(x^2 + y^2 + z^2)/6
%--------------------------------------------------------------------------
f = @(X) ones(size(X,1),1);
gD = @(X) -(X(:,1).^2 + X(:,2).^2 + X(:,3).^2)/6;

% %--------------------------------------------------------------------------
% %   -\delta u = \sin(x)\cos(y)z(1+alpha^2\pi)(1+\alpha^2\pi^2)  ,
% %           u = \sin(x)\cos(y)z(1+alpha^2\pi)
% %--------------------------------------------------------------------------
% alpha = 2;
% f = @(X) sin(X(:,1)).*cos(alpha*pi*X(:,2)).*X(:,3)*(1+alpha^2*pi^2);
% gD = @(X) sin(X(:,1)).*cos(alpha*pi*X(:,2)).*X(:,3);


G = computeGeometry(G);
G = mrstGridWithFullMappings(G);
G = computeVEMGeometry(G,f);

boundaryFaces = (1:G.faces.num)';
boundaryFaces = boundaryFaces( G.faces.centroids(:,1) == 0 | ...
                               G.faces.centroids(:,1) == 1 | ...
                               G.faces.centroids(:,2) == 0 | ...
                               G.faces.centroids(:,2) == 1 | ...
                               G.faces.centroids(:,3) == 0 | ...
                               G.faces.centroids(:,3) == 1 );                          

bc = struct('bcFunc', {{gD}}, 'bcFaces', {{boundaryFaces}}, 'bcType', {{'dir'}});

[A,b] = VEM3D_glob(G,f,bc);

U = A\b;

nodeValues = full(U(1:G.nodes.num));
edgeMidValues = full(U(G.nodes.num + 1:G.nodes.num + G.edges.num));
faceAvg = full(U((G.nodes.num + G.edges.num + 1):(G.nodes.num + G.edges.num + G.faces.num)));
cellAvg = full(U(G.nodes.num + G.edges.num + G.faces.num + 1:end));

figure();
plotFaces(G, 1:G.faces.num, faceAvg);
colorbar;
view(3);

IF = polygonInt3D(G,1:G.faces.num,gD);
IC = polyhedronInt(G,1:G.cells.num,gD);

u = [gD([G.nodes.coords; G.edges.centroids]); IF./G.faces.areas; IC./G.cells.volumes];

err = abs((U - u));
fprintf('Error: %d\n', norm(err, inf));
figure()
plot(err);  
% % plot(nodeValues)
% 
% %   Implement: efficient rule for faceIntegrals.
% %              change bc to give avg values.