function [Sl, bl, dofVec, hK] = VEM3D_loc_v2(G, K)
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
[faceNodes, X] = nodeData(G,K);
                            %   Cell faces, face midpoints and normals.
[faces, Fc, faceNormals] = faceData3D(G,K);
                            %   Cell edges, edge midpoints and normals.
                            %   Baricenter of K.

edgeNum = mcolon(G.faces.edgePos(faces),G.faces.edgePos(faces+1)-1);
                            
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
               X(:,1).^2          , ...   %   (2,0,0)
               X(:,1).*X(:,2) , ...   %   (1,1,0)
               X(:,2).^2 , ...   %   (0,2,0) 
               ];

grad_m3D = @(X) ...
[ones(size(X,1),1)   , zeros(size(X,1),1), zeros(size(X,1),1) ;...    %   (1,0,0)
zeros(size(X,1),1)   , ones(size(X,1),1) , zeros(size(X,1),1) ;...    %   (0,1,0)
zeros(size(X,1),1)   , zeros(size(X,1),1), ones(size(X,1),1)  ;...    %   (0,0,1)
X(:,1)*2             , zeros(size(X,1),1), zeros(size(X,1),1) ;...    %   (2,0,0)
X(:,2)               , X(:,1)            , zeros(size(X,1),1) ;...    %   (1,1,0)
X(:,3)               , zeros(size(X,1),1), X(:,1)             ;...    %   (1,0,1)
zeros(size(X,1),1)   , X(:,2)*2          , zeros(size(X,1),1) ;...    %   (0,2,0)
zeros(size(X,1),1)   , X(:,3)            , X(:,2)             ;...    %  (0,1,1)
zeros(size(X,1),1)   , zeros(size(X,1),1), X(:,3)*2];    %   (0,0,2)

Kc = G.cells.centroids(K,:);

hK = G.cells.diameters(K);     %   Element diameter.
vol = G.cells.volumes(K);   %   Element volume

nV = size(X,1);             %   Number of nodes.
nF = size(faces,1);         %   Number of faces.

k = 2;                      %   Method order.
nk = (k+1)*(k+2)*(k+3)/6;   %   Dimension of polynom[1,G.faces.faceInt{1}(faces(i),:)]'ial space.
                            %   Local nomber of dofs.
% NK = nV + nE*(k-1) + nF*k*(k-1)/2 + k*(k^2-1)/6;

%%  BUILD FACE MATRICES

Fc = G.faces.centroids(faces,:);
hF = G.faces.diameters(faces,:);

grad_m = @(X) [...
    ones(size(X,1),1)    , zeros(size(X,1),1), ...
    zeros(size(X,1),1)   , ones(size(X,1),1) , ...
    X(:,1)*2             , zeros(size(X,1),1), ...
    X(:,2)               , X(:,1)            , ...
    zeros(size(X,1),1)   , X(:,2)*2          ];

int_m = @(X)        [X(:,1).^2/2  , ...
                     X(:,1).*X(:,2)   , ...
                     X(:,1).^3/3         , ...
                     X(:,1).^2.*X(:,2)/2 , ...
                     X(:,1).*X(:,2).^2      , ...
                     ];

for i = 1:nF
    
    edgeNum = G.faces.edgePos(faces(i)):G.faces.edgePos(faces(i)+1)-1;
    edges = G.faces.edges(edgeNum);
    edgeNormals = G.faces.edgeNormals(edgeNum,:);
    
    nodeNum = mcolon(G.edges.nodePos(edges),G.edges.nodePos(edges+ 1)-1);
    faceNodes = G.edges.nodes(nodeNum);
    faceNodes = reshape(faceNodes,2,[])';
    nN = size(faceNodes,1);
    
    
    NF = 2*nN+1;
    faceNodes(G.faces.edgeSign(edgeNum) == -1,:) = ...
                             faceNodes(G.faces.edgeSign(edgeNum) == -1,2:-1:1);
    faceNodes = reshape(faceNodes,[],1);
    vec1 = (G.nodes.coords(faceNodes(nN+1),:) - G.nodes.coords(faceNodes(1),:))';
                            %   NB: faceNormal must be in correct
                            %   direction, as done by faceData.
    vec2 = cross(G.faces.normals(i,:)',vec1);
    T = [vec1,vec2];
    XF = [G.nodes.coords(faceNodes,:);G.edges.centroids(edges,:)];
    XF = bsxfun(@rdivide,XF-repmat(Fc(i,:),size(XF,1),1), ...
                                 hF(i).*ones(size(XF,1),1));
    XF = XF*T;
    edgeNormals = edgeNormals*T;
    edgeNormals = bsxfun(@times, edgeNormals,G.edges.lengths(edges));
    I = bsxfun(@times,[(int_m(XF(1:nN,:)) + int_m(XF(nN+1:2*nN,:)))/6 ... 
           + int_m(XF(2*nN+1:end,:))*2/3], edgeNormals(:,1));
    I = sum(I, 1).*hF(i);
                           %   Build matrix Df
    DF = [ones(2*nN,1), m(XF([1:nN, 2*nN+1:3*nN],:)); [1,I]];
    
                            %   Build matrix Bf
    BF = zeros(6, 2*nN+1);
    BF(1,NF) = 1;
    vals = grad_m(XF).*repmat(edgeNormals,3,5);
                                %   Fix this...
    vals = [sum(vals(:,1:2),2), sum(vals(:,3:4),2),sum(vals(:,5:6),2), ...
            sum(vals(:,7:8),2), sum(vals(:,9:10),2)];

    I = [(vals(1:nN,:) + vals([2*nN,nN+1:2*nN-1],:))/6; ...
                         vals(2*nN+1:end,:)*2/3]./hF(i);
        
    BF(2:6,1:NF-1) = I';
    BF([4,6],NF) = -2*G.faces.areas(faces(i))/G.faces.diameters(faces(i)).^2;
    MF = BF*DF;
    PNFstar = MF\BF;
    PNF = DF*PNFstar;
    f = @(X) X(:,1).^2 + X(:,3).^2 - X(:,2)/8;
%     fv = [f(X); polygonInt(T*X(1:nN,:),f)/G.faces.areas(i)];
    
    XKF = [G.nodes.coords(faceNodes,:);G.edges.centroids(edges,:)];
    XKF = bsxfun(@rdivide,XKF - repmat(Kc,size(XKF,1),1), ...
                                 hK.*ones(size(XF,1),1));
    
    vals = sum(grad_m3D(XKF).*repmat(faceNormals(i,:),27*nN,1),2);
    vals = reshape(vals, 3*nN, 9)
    
    
end



end
