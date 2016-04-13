%--------------------------------------------------------------------------
%   Minimal working example of VEM in three dimensions.
%--------------------------------------------------------------ØSK-2016021-

clc; clear; close all;

addpath('../')              %  Extra grid mappings.
            
                            %   Voronoi grid generation. Replace path with
                            %   your own local version.
addpath('/home/strene/Documents/master/coop/pebiGridding/voronoi3D')
k = 2;
                            %   Generate grid.
G = voronoiCube(50,[1,1,1]);

                            %   Set source term
f = @(X) ones(size(X,1),1);

                            %   Compute VEM geometry.
G = computeVEM3DGeometry(G);

                            %   Dirichlet boundary condition.
                            %   This is also the exact solution of this
                            %   example.
gD = @(X) -(X(:,1).^2 + X(:,2).^2 + X(:,3).^2)/6;

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

                            %   Plot face avarages, with ball of radius r,
                            %   center (0,0,0), carved out.
Kc = G.cells.centroids;
cells = 1:G.cells.num;
r = .7;
cells = cells(Kc(:,1).^2 + Kc(:,2).^2 + Kc(:,3).^2 > r^2);
faceNum = mcolon(G.cells.facePos(cells),G.cells.facePos(cells+1)-1);
faces = G.cells.faces(faceNum);

figure();
plotFaces(G,faces,sol.faceMoments(faces));
colorbar;
view(3);
axis equal;