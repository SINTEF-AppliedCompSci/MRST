%--------------------------------------------------------------------------
%   Convergence test.
%--------------------------------------------------------------ØSK-2016021-

clc; clear; close all;

addpath('../')              %  Extra grid mappings.
            
                            %   Voronoi grid generation. Replace path with
                            %   your own local version.
addpath('/home/strene/Documents/master/coop/pebiGridding/voronoi3D')


% f = @(X) 100*X(:,1).*pi^2.*cos(pi*X(:,2)) - exp(X(:,3));
% 
%                             %   Dirichlet boundary condition.
%                             %   This is also the exact solution of this
%                             %   example.
% gD = @(X) 100*X(:,1).*cos(pi*X(:,2)) + exp(X(:,3));

% %--------------------------------------------------------------------------
% %   -\delta u = 0,
% %           u = 1/(2\pi||x-C||)
% %--------------------------------------------------------------------------
% f = @(X) zeros(size(X,1),1);
% C = -[.2,.2,.2];
% gD = @(X) -1./(2*pi*sqrt(sum((X-repmat(C,size(X,1),1)).^2,2)));

f = @(X) -exp(X(:,1));
gD = @(X) exp(X(:,1));

nVec = [200 400 800 1600];
% nVec = [40, 160, 640];
% nVec = [10, 20, 30, 40];
nn = numel(nVec);
hVec = zeros(nn, 1);
errVec = zeros(nn,1);

for i = 1:nn
    
n = nVec(i);
                            %   Generate grid.
G = voronoiCube(n,[1,1,1]);
n = ceil(n^(1/3));
% G = cartGrid([n,n,n], [1,1,1]);

                            %   Compute VEM geometry.
G = computeVEMGeometry(G,f);

                            %   Define boundary faces.
boundaryFaces = (1:G.faces.num)';                           
boundaryFaces = boundaryFaces( G.faces.neighbors(:,1) == 0 | ...
                               G.faces.neighbors(:,2) == 0       );

                            %   Set boundary condition struct.
bc = struct('bcFunc' , {{gD}}           , ...
            'bcFaces', {{boundaryFaces}}, ...
            'bcType' , {{'dir'}}              );

                            %   Solve using second order VEM.
sol = VEM3D(G,f,bc,2);
U = [sol.nodeValues; sol.edgeValues; sol.faceMoments; sol.cellMoments];

IF = polygonInt3D(G,1:G.faces.num,gD, 7);
IC = polyhedronInt(G,1:G.cells.num,gD, 7);
u = [gD([G.nodes.coords; G.edges.centroids]);
     IF./G.faces.areas ; IC./G.cells.volumes ];

nK = G.cells.num;
h = sum(G.cells.diameters)/nK;
hVec(i) = h;

errVec(i) = sqrt(h^3)*norm(U-u, 2);

end

loglog(hVec, errVec, '-sq', hVec, .8*errVec(nn)*(hVec/hVec(nn)).^3, '-.r', 'LineWidth', 1.5);

slope = polyfit(log(hVec), log(errVec), 1);    
fprintf('Slope: %f', slope(1));