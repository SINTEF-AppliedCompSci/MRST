clc; clear; close all;

addpath('/home/strene/Documents/master/coop/pebiGridding/voronoi3D')

n = 10;
gridLim = [1,1,1];

% G = cartGrid([n,n,n],gridLim);
% G = tetrahedronCube([n,n,n], gridLim, 1);
G = voroniCube(n^3,gridLim);


% G = computeGeometry(G);
% 
% remCells = (1:G.cells.num)';
% remCells = remCells(G.cells.centroids(:,3) > 1 - ...
%                     (G.cells.centroids(:,1) + G.cells.centroids(:,2))   );
% G = removeCells(G,remCells);
% 
% R1 = @(theta) [1 0 0; 0 cos(theta) sin(theta); 0 -sin(theta) cos(theta)];   
% R2 = @(phi)   [cos(phi) 0 -sin(phi); 0 1 0; sin(phi) 0 cos(phi)];
% R3 = @(psi)   [cos(psi) sin(psi) 0; -sin(psi) cos(psi) 0; 0 0 1];  
% 
% X = G.nodes.coords;
% 
% % X(:,2) = X(:,2) + 1/10*X(:,3).^2;
% 
% theta = pi/6; phi = pi/3; psi = pi/18;
% X = X*R1(theta)*R2(phi)*R3(psi);
% 
% G.nodes.coords = X;



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

% for i = 1:G.faces.num
%     edgeNum = G.faces.edgePos(i):G.faces.edgePos(i+1)-1;
%     edges = G.faces.edges(edgeNum);
%     nodeNum = mcolon(G.edges.nodePos(edges),G.edges.nodePos(edges+1)-1);
%     nodes = G.edges.nodes(nodeNum);
%     X = G.nodes.coords(nodes,:);
%     nX = size(X,1);
%     for j = 1:nX/2
%         plot3([X(j,1); X(j+1,1)], ...
%               [X(j,2); X(j+1,2)], ...
%               [X(j,3); X(j+1,3)])
% %         plot3([X(j,1); X(mod(j,nX) + 1, 1)], ...
% %               [X(j,2); X(mod(j,nX) + 1, 2)], ...
% %               [X(j,3); X(mod(j,nX) + 1, 3)]);
%         hold on
%         pause
%     end
%     hold off
% end
    

boundaryFaces = (1:G.faces.num)';
% boundaryFaces = boundaryFaces( G.faces.centroids(:,1) == 0          | ...
%                                G.faces.centroids(:,1) == gridLim(1) | ...
%                                G.faces.centroids(:,2) == 0          | ...
%                                G.faces.centroids(:,2) == gridLim(2) | ...
%                                G.faces.centroids(:,3) == 0          | ...
%                                G.faces.centroids(:,3) == gridLim(3) );
                           
boundaryFaces = boundaryFaces( G.faces.neighbors(:,1) == 0 | ...
                               G.faces.neighbors(:,2) == 0 );

bc = struct('bcFunc', {{gD}}, 'bcFaces', {{boundaryFaces}}, 'bcType', {{'dir'}});

[A,b] = VEM3D_glob(G,f,bc);

U = A\b;

clear A b

nodeValues = full(U(1:G.nodes.num));
edgeMidValues = full(U(G.nodes.num + 1:G.nodes.num + G.edges.num));
faceAvg = full(U((G.nodes.num + G.edges.num + 1):(G.nodes.num + G.edges.num + G.faces.num)));
cellAvg = full(U(G.nodes.num + G.edges.num + G.faces.num + 1:end));


Kc = G.cells.centroids;
cells = 1:G.cells.num;
r = .7;
cells = cells(Kc(:,1).^2 + Kc(:,2).^2 + Kc(:,3).^2 > r^2);
faceNum = mcolon(G.cells.facePos(cells),G.cells.facePos(cells+1)-1);
faces = G.cells.faces(faceNum);
% faces = 1:G.faces.num;

figure();
plotFaces(G,faces,faceAvg(faces));

colorbar;
view(3);
axis equal;

IF = polygonInt3D(G,1:G.faces.num,gD);
IC = polyhedronInt(G,1:G.cells.num,gD);

u = [gD([G.nodes.coords; G.edges.centroids]); IF./G.faces.areas; IC./G.cells.volumes];

err = abs((U - u));

h = sum(G.cells.diameters)/G.cells.num;

fprintf('Error: %d\n', h^(3/2)*norm(err, 2));
figure()
plot(err);  
% % plot(nodeValues)
% 
% %   Implement: efficient rule for faceIntegrals.
% %              change bc to give avg values.




% vols1 = zeros(G.cells.num,1);
% vols2 = zeros(G.cells.num,1);
% 
% 
% V  = [-1.0, -1.0/sqrt(3.0), -1.0/sqrt(6.0); ...
%        0.0,  2.0/sqrt(3.0), -1.0/sqrt(6.0); ...
%        1.0, -1.0/sqrt(3.0), -1.0/sqrt(6.0); ...
%        0.0,  0.0          ,  3.0/sqrt(6.0)];   
% Vdiff = V(1:end-1,:) - V(2:end,:);
%    
% vol = sqrt(sqrt(8.0)/3.0);
% 
% for i = 1:G.cells.num
%     nodeNum = G.cells.nodePos(i):G.cells.nodePos(i+1)-1;
%     nodes = G.cells.nodes(nodeNum);
%     X = G.nodes.coords(nodes,:);
%     tri = delaunay(X);
%     nTri = size(tri,1);
%     for j = 1:nTri
%         Xdiff = X(tri(j,1:3),:) - X(tri(j,2:4),:);
%         A1 = Vdiff\Xdiff;
%         b = X(tri(j,1),:) - V(1,:)*A1; 
%         Xhat = V*A1 + repmat(b,size(V,1),1);
%         X(tri(j,:),:) - Xhat;
%         
%         A2 = X(tri(j,2:4),:) - repmat(X(tri(j,1),:),3,1);
%         
%         vols1(i) = vols1(i) + abs(det(A1))*vol;
%         vols2(i) = vols2(i) + abs(det(A2))/6;
%         
%     end
% end
% G.cells.volumes - vols1;
% G.cells.volumes - vols2;
% 
% % G.cells.volumes = vols2;
