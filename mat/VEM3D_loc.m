function [Sl, bl, dofVec, hK] = VEM3D_loc(G, K)
%--------------------------------------------------------------------------
%   Generates local stffness matrix for the virtual element method  for
%   cell K of grid G for diffusion problem:
%
%           -\delta u = f, x \in \Omega
%                   u = g, x \in \partial \Omega
%
%   Input:
%
%   G:      2D MRST grid. Cells can be any kind of polygn. the function
%           assumes the following functions has been called for the grid:
%
%           G = mrstGridWithFullMappings(G);
%           G = computeGeometry(G);
%           G = globalEdgeData(G);
%
%   K:      Cell number in grid G, i.e. G.cells(K).
%
%   f:      Source term.
%
%   Output:
%
%   Sl:     Local stiffness matrix for cell K. dim(Sl) = NK x NK, where
%           NK = n*k + 0.5*k*(k-1), n is the number of vertices of K,
%           and k = 2 is the order of the method.
%
%   bl:     Local load vector for cell K. dim(bl) = 1 x NK.
%   
%   dofVec: Map from local to global dofs. S(dofVec, dofVec) = Sl, where
%           S is the global stiffness matrix.
% 
%--------------------------------------------------------------------------

%%  CELL DATA                                                            %%


                            %   Cell nodes and node coordinates.
[nodes, X] = nodeData(G,K);
                            %   Cell faces, face midpoints and normals.
[faces, Fc, faceNormals] = faceData(G,K);
                            %   Cell edges, edge midpoints and normals.
                            %   Baricenter of K.
% nodeNum = mcolon(G.edges.nodePos(edges),G.edges.nodePos(edges+1)-ones(size(edges,1),1))
% nodes = G.edges.nodes(nodeNum)
% plot3(G.nodes.coords(nodes,1), G.nodes.coords(nodes,2), G.nodes.coords(nodes,3), '*')

                            %   m =      @(X) [(X(:,1)-xK)/hK               , ...   %   (1,0,0)
%                (X(:,2)-yK)/hK               , ...   %   (0,1,0)
%                (X(:,3)-zK)/hK               , ...   %   (0,0,1)
%                (X(:,1)-xK).^2/hK^2          , ...   %   (2,0,0)
%                (X(:,1)-xK)*(X(:,2)-yK)/hK^2 , ...   %   (1,1,0)
%                (X(:,1)-xK)*(X(:,3)-zK)/hK^2 , ...   %   (1,0,1)
%                (X(:,2)-xK).^2/hK^2          , ...   %   (0,2,0) 
%                (X(:,2)-yK)*(X(:,3)-zK)/hK^2 , ...   %   (0,1,1)
%                (X(:,3)-zK).^2/hK^2          ];      %   (0,0,2)icenter of K.

m =      @(X) [X(:,1)               , ...   %   (1,0,0)
               X(:,2)               , ...   %   (0,1,0)
               X(:,3)               , ...   %   (0,0,1)
               X(:,1).^2          , ...   %   (2,0,0)
               X(:,1).*X(:,2) , ...   %   (1,1,0)
               X(:,1).*X(:,3) , ...   %   (1,0,1)
               X(:,2).^2 , ...   %   (0,2,0) 
               X(:,2).*X(:,3) , ...   %   (0,1,1)
               X(:,3).^2          ];      %   (0,0,2)

Kc = G.cells.centroids(K,:);

hK = G.cells.diameters(K);     %   Element diameter.
vol = G.cells.volumes(K);   %   Element volume

nV = size(X,1);             %   Number of nodes.
nF = size(faces,1);         %   Number of faces.

k = 2;                      %   Method order.
nk = (k+1)*(k+2)*(k+3)/6;   %   Dimension of polynomial space.
                            %   Local nomber of dofs.
%NK = nV + nE*(k-1) + nF*k*(k-1)/2 + k*(k^2-1);

%%  BUILD FACE MATRICES

diffVec = diff(G.faces.nodePos);
Fc = G.faces.centroids(faces,:);
hF = G.faces.diameters(faces,:);

for i = 1:nF
    faceNum = G.faces.edgePos(faces(i)):G.faces.edgePos(faces(i)+1)-1;
    faceEdges = G.faces.edges(faceNum);
    nodeNum = G.faces.nodePos(faces(i)):G.faces.nodePos(faces(i) + 1)-1;
    faceNodes  = G.faces.nodes(nodeNum);
    Xf = [G.nodes.coords(faceNodes,:);G.edges.centroids(faceEdges,:)];
    Xf = bsxfun(@rdivide,Xf-repmat(Fc(i,:),size(Xf,1),1), ...
                                 hF(i).*ones(size(Xf,1),1));
                            %   Build matrix Df
    Df = [ones(size(Xf,1),1), m(Xf); ...
         [1,G.faces.faceInt{1}(faces(i),:)]];
                            %   Build matrix Bf
    vec = zeros(9,1);
    vec([4,7,9]) = -2*G.faces.areas(i)/G.faces.diameters(i).^2;
    Bf = [zeros(1,2*diffVec(i)), 1 ; ...
         G.faces.faceInt{2}([faceNum, faceNum + 2*diffVec(i)],:)', vec];
    Gf = Bf*Df;
end