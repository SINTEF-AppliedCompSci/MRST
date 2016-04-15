clc; clear; close all;

addpath('../')
addpath('/home/strene/Documents/master/coop/pebiGridding/voronoi3D')

n = 8;
gridLim = [1,1,1];

% G = cartGrid([n,n,n],gridLim);
% 
% X = G.nodes.coords;
% X = bsxfun(@minus, X, [1,1,1]);

% G.nodes.coords = X;
% G = tetrahedronCube([n,n,n], gridLim, 1);
G = voronoiCube(100  ,gridLim);


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

%--------------------------------------------------------------------------
%   -\delta u = 0,
%           u = 1/(2\pi||x-C||)
%--------------------------------------------------------------------------
f = @(X) -2*ones(size(X,1),1);
gD = @(X) X(:,1).^2;
k = 2;

G = computeVEM3DGeometry(G);
                           
boundaryFaces = find( G.faces.neighbors(:,1) == 0 | ...
                      G.faces.neighbors(:,2) == 0 );

bc = VEM3D_addBC([], boundaryFaces, 'pressure', gD);

[sol, G] = VEM3D(G, f, k, bc, 'cellProjectors', true);
U = [sol.nodeValues; sol.edgeValues; sol.faceMoments; sol.cellMoments];

Kc = G.cells.centroids;
cells = 1:G.cells.num;
r = .7; c = [1,0,0];
cells = cells(sum(bsxfun(@minus, Kc, c).^2,2) > r^2);
faceNum = mcolon(G.cells.facePos(cells),G.cells.facePos(cells+1)-1);
faces = G.cells.faces(faceNum);

if k == 2

figure();
plotFaces(G,faces,sol.faceMoments(faces));
colorbar;
view(3);
axis equal;

IF = polygonInt3D(G,1:G.faces.num,gD, 7);
IC = polyhedronInt(G,1:G.cells.num,gD, 7);

u = [gD([G.nodes.coords; G.edges.centroids]); IF./G.faces.areas; IC./G.cells.volumes];
err = abs((U - u));
elseif k == 1
    u = gD(G.nodes.coords);
    err = abs(U-u);
end

l2Err = l2Error3D(G, sol, gD ,k);

fprintf('2-norm error: %d\n', norm(err, 2));
fprintf('L^2-norm error: %d\n', sqrt(sum(l2Err.^2)));
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
