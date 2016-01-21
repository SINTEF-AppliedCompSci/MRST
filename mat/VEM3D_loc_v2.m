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
[nodes, X] = nodeData(G,K);
nN = size(nodes,1);


                            %   Cell faces, face midpoints and normals.
[faces, Fc, faceNormals] = faceData3D(G,K);
                            %   Cell edges, edge midpoints and normals.
                            %   Baricenter of K.

                            % nodeNum = mcolon(G.edges.nodePos(edges),G.edges.nodePos(edges+1)-ones(size(edges,1),1))
% nodes = G.edges.nodes(nodeNum)
% plot3(G.nodes.coords(nodes,1), G.nodes.coords(nodes,2), G.nodes.coords(nodes,3), '*')

m3D =      @(X) [ones(size(X,1),1) , ...
                X(:,1)              , ...   %   (1,0,0)
               X(:,2)               , ...   %   (0,1,0)
               X(:,3)               , ...   %   (0,0,1)
               X(:,1).^2          , ...   %   (2,0,0)
               X(:,1).*X(:,2) , ...   %   (1,1,0)
               X(:,1).*X(:,3), ...   %   (1,0,1)
               X(:,2).^2, ...   %   (0,2,0) 
               X(:,2).*X(:,3), ...   %   (0,1,1)
               X(:,3).^2];      %   (0,0,2)icenter of K.

m2D =      @(X) [ones(size(X,1),1), ...
               X(:,1)               , ...   %   (1,0,0)
               X(:,2)               , ...   %   (0,1,0)
               X(:,1).^2          , ...   %   (2,0,0)
               X(:,1).*X(:,2) , ...   %   (1,1,0)
               X(:,2).^2 , ...   %   (0,2,0) 
               ];

int_m3D = @(X)        [X(:,1).^2/2  , ...
                     X(:,1).*X(:,2)   , ...
                     X(:,1).*X(:,3)   , ...
                     X(:,1).^3/3         , ...
                     X(:,1).^2.*X(:,2)/2 , ...
                     X(:,1).^2.*X(:,3)/2  , ...
                     X(:,1).*X(:,2).^2      , ...
                     X(:,1).*X(:,2).*X(:,3) , ...
                     X(:,1).*X(:,3).^2      ];
           
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

% grad_m3D = @(X) ...
% [X(:,1)*2            , zeros(size(X,1),1), zeros(size(X,1),1) ;...    %   (2,0,0)
% X(:,2)               , X(:,1)            , zeros(size(X,1),1) ;...    %   (1,1,0)
% X(:,3)               , zeros(size(X,1),1), X(:,1)             ;...    %   (1,0,1)
% zeros(size(X,1),1)   , X(:,2)*2          , zeros(size(X,1),1) ;...    %   (0,2,0)
% zeros(size(X,1),1)   , X(:,3)            , X(:,2)             ;...    %  (0,1,1)
% zeros(size(X,1),1)   , zeros(size(X,1),1), X(:,3)*2];    %   (0,0,2)


edgeNum = G.cells.edgePos(K):G.cells.edgePos(K+1)-1;
edges = G.cells.edges(edgeNum);
nE = size(edges,1);

X = [G.nodes.coords(nodes,:);
     G.edges.centroids(edges,:);
     G.faces.centroids(faces,:)];

Kc = G.cells.centroids(K,:);

hK = G.cells.diameters(K);     %   Element diameter.

X = (X-repmat(Kc,size(X,1),1))./hK;

vol = G.cells.volumes(K);   %   Element volume

nV = size(X,1);             %   Number of nodes.
nF = size(faces,1);         %   Number of faces.

k = 2;                      %   Method order.
nk = (k+1)*(k+2)*(k+3)/6;   %   Dimension of polynom[1,G.faces.faceInt{1}(faces(i),:)]'ial space.
                            %   Local nomber of dofs.
% NK = nV + nE*(k-1) + nF*k*(k-1)/2 + k*(k^2-1)/6;

NK = nN + nE*(k-1) + nF*k*(k-1)/2 + k*(k^2-1)/6;

%%  BUILD FACE MATRICES

Fc = G.faces.centroids(faces,:);
hF = G.faces.diameters(faces,:);
faceAreas = G.faces.areas(faces);

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

B = zeros(nk, NK);
D = zeros(NK,nk);
D(1:NK-1,:) = m3D(X(1:NK-1,:));
                 
for i = 1:nF
    
    edgeNum = G.faces.edgePos(faces(i)):G.faces.edgePos(faces(i)+1)-1;
    faceEdges = G.faces.edges(edgeNum);
    edgeNormals = G.faces.edgeNormals(edgeNum,:);
    
    nodeNum = mcolon(G.edges.nodePos(faceEdges),G.edges.nodePos(faceEdges+ 1)-1);
    faceNodes = G.edges.nodes(nodeNum);
    faceNodes = reshape(faceNodes,2,[])';
    nNF = size(faceNodes,1);
    NF = 2*nNF+1;
    faceNodes(G.faces.edgeSign(edgeNum) == -1,:) = ...
        faceNodes(G.faces.edgeSign(edgeNum) == -1,2:-1:1);
    faceNodes = reshape(faceNodes,[],1);
    vec1 = (G.nodes.coords(faceNodes(nNF+1),:) - G.nodes.coords(faceNodes(1),:))';
    vec1 = vec1/norm(vec1,2);
                            %   NB: faceNormal must be in correct
                            %   direction, as done by faceData.
    vec2 = cross(faceNormals(i,:)',vec1);
    vec2 = vec2/norm(vec2,2);
    T = [vec1,vec2];
    XFmon = [G.nodes.coords(faceNodes,:);G.edges.centroids(faceEdges,:)];
    
    XF = (XFmon-repmat(Fc(i,:),size(XFmon,1),1))*T;
    XFmon = bsxfun(@rdivide,XFmon-repmat(Fc(i,:),size(XFmon,1),1), ...
                                 hF(i).*ones(size(XFmon,1),1));

    XFmon = XFmon*T;
    edgeNormals = edgeNormals*T;
    edgeNormals = bsxfun(@times, edgeNormals,G.edges.lengths(faceEdges));
    I = bsxfun(@times,[(int_m(XFmon(1:nNF,:)) + int_m(XFmon(nNF+1:2*nNF,:)))/6 ... 
           + int_m(XFmon(2*nNF+1:end,:))*2/3], edgeNormals(:,1));
    I = sum(I, 1).*hF(i);
    ID = I;    
                           %   Build matrix Df
    DF = [m2D(XFmon([1:nNF, 2*nNF+1:3*nNF],:)); [1,I/faceAreas(i)]];
    
                            %   Build matrix Bf
    BF = zeros(6, 2*nNF+1);
    BF(1,NF) = 1;
    vals = grad_m(XFmon).*repmat(edgeNormals,3,5);
                                %   Fix this...
    vals = [sum(vals(:,1:2),2), sum(vals(:,3:4),2),sum(vals(:,5:6),2), ...
            sum(vals(:,7:8),2), sum(vals(:,9:10),2)];

    I = [(vals(1:nNF,:) + vals([2*nNF,nNF+1:2*nNF-1],:))/6; ...
                          vals(2*nNF+1:end,:)*2/3]./hF(i);
        
    BF(2:6,1:NF-1) = I';
    BF([4,6],NF) = -2*faceAreas(i)/hF(i)^2;
% 
%     BFT = BF; BFT(:, [5,8]) = BFT(:, [8,5]);
%     DFT = DF; DFT([5,8],:) = DFT([8,5],:);
%     
%     DF = DFT; BF  =BFT;
    
    MF = BF*DF;
%     DFT = DF; DFT([5,8],:) = DFT([8,5],:);
    PNFstar = MF\BF;
    PNF = DF*PNFstar;
%     f = @(X) X(:,1).^2  + X(:,2)/8;
%     ii= [1:nNF, 2*nNF+1:3*nNF];
%     ii = [1:nNF, 3*nNF, 2*nNF+1:3*nNF-1];
%     fv = [f(XF(ii,:)); polygonInt(XF(1:nNF,:),f)/faceAreas(i)];
%     fv - [m2D(XFmon(ii,:)); 1, ID/faceAreas(i)] *PNFstar*fv
%     fv - PNF*fv
    
    XKF = [G.nodes.coords(faceNodes,:);G.edges.centroids(faceEdges,:)];
    XKF = bsxfun(@rdivide,XKF - repmat(Kc,size(XKF,1),1), ...
                                 hK.*ones(size(XFmon,1),1))*T;
    
    I = polygonFaceInt(XF(1:nNF,:),hK, Kc, hF(i), Fc(i,:), T, faceNormals(i,:)./faceAreas(i), ...
                              grad_m3D, m2D, int_m3D, PNFstar);
    
    [~, iiN] =  ismember(faceNodes(1:nNF), nodes);
    [~, iiE] = ismember(faceEdges,edges);
   
    dofVec = [iiN', iiE' + nN, i + nN + nE];
    
     B(2:10,dofVec) = B(2:10,dofVec) + I;
       
                            %   NB: Divergence rule might improve speed.
    
    end

B(1,NK) = 1;
B([5,8,10],NK) = -2*G.faces.areas(faces(i))/G.faces.diameters(faces(i)).^2;


D = [m3D(X); G.cells.monomialIntegrals(K,:)./vol];
M = B*D
end
