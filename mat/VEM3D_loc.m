function [Sl, bl, dofVec] = VEM3D_loc_v2(G, f, K)
%--------------------------------------------------------------------------
%
%   Awesome function.
%
%--------------------------------------------------------------------------

%%  CELL DATA                                                            %%

                            %   Node data for cell K.
nodeNum = G.cells.nodePos(K):G.cells.nodePos(K+1)-1;
nodes   = G.cells.nodes(nodeNum);
X       = G.nodes.coords(nodes,:);
nN      = size(nodes,1);

                            %   Edge data for cell K.
edgeNum = G.cells.edgePos(K):G.cells.edgePos(K+1)-1;
edges   = G.cells.edges(edgeNum);
Ec      = G.edges.centroids(edges,:);
nE      = size(edges,1);

                            %   Face data for cell K.
faceNum     = G.cells.facePos(K):G.cells.facePos(K+1)-1;
faces       = G.cells.faces(faceNum);
Fc          = G.faces.centroids(faces,:);
faceNormals = G.faces.normals(faces,:);
nF          = size(faces,1);

                            %   Cell data for cell K.
Kc  = G.cells.centroids(K,:);
hK  = G.cells.diameters(K);
vol = G.cells.volumes(K);

k  = 2;                     %   Method order.
nk = (k+1)*(k+2)*(k+3)/6;   %   Dimension of polynomial space.
                            %   Local nomber of dofs.
NK = nN + nE*(k-1) + nF*k*(k-1)/2 + k*(k^2-1)/6;

X = (X-repmat(Kc,size(X,1),1))/hK;

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

for i = 1:nF
    
    PNFstarNum = G.faces.PNFstarPos(i):G.faces.PNFstar(i+1)-1;
    PNFstar = G.faces.PNFstar(PNFstarNum,:);
    
    XKF = [G.nodes.coords(faceNodes,:);G.edges.centroids(faceEdges,:)];
    XKF = bsxfun(@rdivide,XKF - repmat(Kc,size(XKF,1),1), ...
                                 hK.*ones(size(XFmon,1),1))*T;
    
    I = polygonFaceInt(XF(1:nNF,:),hK, Kc, hF(i), Fc(i,:), T, ...
        faceNormals(i,:)./faceAreas(i), ...
                              grad_m3D, m2D, int_m3D, PNFstar);
    
    [~, iiN] =  ismember(faceNodes(1:nNF), nodes);
    [~, iiE] = ismember(faceEdges,edges);
   
    dofVec = [iiN', iiE' + nN, i + nN + nE];
    
    B(2:10,dofVec) = B(2:10,dofVec) + I;
    
end

B(1,NK) = 1;
B([5,8,10],NK) = -2*vol/hK.^2;
faceIntNum = G.cells.faceIntPos(K):G.cells.faceIntPos(K+1)-1;
D = [m3D(X); ...
     bsxfun(@rdivide, G.cells.monomialFaceIntegrals(faceIntNum,:), ...
                      G.faces.areas(faces)); ...
     G.cells.monomialCellIntegrals(K,:)./vol];
M = B*D;

PNstar = M\B;
PN = D*PNstar;
% 
% f = @(X) X(:,1).^2 + X(:,3).*X(:,2)/27*60;
% 
% fF = polygonInt3D(G,f);
% fK = polyhedronInt2(G,f);
% 
% X = [G.nodes.coords(nodes,:); G.edges.centroids(edges,:)];
% 
% fv = [f(X); fF(faces)./G.faces.areas(faces); fK(K)/vol];
% 
% fv - PN*fv

Mtilde = [zeros(1,nk); M(2:nk,:)];

Sl = PNstar'*Mtilde*PNstar + hK*(eye(NK)-PN)'*(eye(NK)-PN);

%%  LOCAL LOAD TERM.                                                     %%

                            %   Matrix of integrals over K of all
                            %   combinations of linear monomials.
I = G.cells.monomialCellIntegrals(K,:); 

H = [I(1:4)        ; ...
     I([2,5,6,7])  ; ...
     I([3,6,8,9])  ; ...
     I([4,7,9,10])];             
                            %   \Pi^\Nabla in the monomial basis
                            %   \mathcal{M}_1.
PNstar = M(1:4,1:4)\B(1:4,:);
                            %   Local load term.
                            
fChi = [f([G.nodes.coords(nodes,:); G.edges.centroids(edges,:)]); ...
        G.faces.fFaceIntegrals(faces)./G.faces.areas(faces)                 ; ...
        G.cells.fCellIntegrals(K)/vol];

% fChi - PN*fChi

bl = PNstar'*H*PNstar*fChi;


dofVec = [nodes', edges' + G.nodes.num, ...
          faces' + G.nodes.num + G.edges.num, ...
          K + G.nodes.num + G.edges.num + G.faces.num];
      
end