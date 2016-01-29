function [Sl, bl, dofVec] = VEM3D_loc(G, f, K)
%--------------------------------------------------------------------------
%
%   Awesome function.
%
%--------------------------------------------------------------------------

%%  CELL DATA                                                            %%

                            %   Node data for cell K.
nodeNum = G.cells.nodePos(K):G.cells.nodePos(K+1)-1;
nodes   = G.cells.nodes(nodeNum);
if size(nodes,1) == 1;
    nodes = nodes';
end
X       = G.nodes.coords(nodes,:);
monNodeNum = G.cells.monomialNodeValsPos(K):G.cells.monomialNodeValsPos(K+1)-1;
monomialNodeVals = G.cells.monomialNodeVals(monNodeNum,:);
nN      = size(nodes,1);

                            %   Edge data for cell K.
edgeNum = G.cells.edgePos(K):G.cells.edgePos(K+1)-1;
edges   = G.cells.edges(edgeNum);
if size(edges,1) == 1;
    edges = edges';
end
Ec      = G.edges.centroids(edges,:);
monEdgeNum = G.cells.monomialEdgeValsPos(K):G.cells.monomialEdgeValsPos(K+1)-1;
monomialEdgeVals = G.cells.monomialEdgeVals(monEdgeNum,:);
nE      = size(edges,1);

                            %   Face data for cell K.
faceNum     = G.cells.facePos(K):G.cells.facePos(K+1)-1;
faces       = G.cells.faces(faceNum);
if size(faces,1) == 1;
    faces = faces';
end
faceAreas   = G.faces.areas(faces);
nF          = size(faces,1);
faceIntNum = G.cells.faceIntPos(K):G.cells.faceIntPos(K+1)-1;
monomialFaceInt     = bsxfun(@rdivide, ...
                     G.cells.monomialFaceIntegrals(faceIntNum,:), ...
                     faceAreas);
mIint = G.cells.monomialFaceIntegrals(faceIntNum,:)
fFaceIntegrals = G.faces.fFaceIntegrals(faces);

                 
                            %   Cell data for cell K.
Kc  = G.cells.centroids(K,:);
hK  = G.cells.diameters(K);
vol = G.cells.volumes(K);
monomialCellInt = G.cells.monomialCellIntegrals(K,:);
fCellIntegral = G.cells.fCellIntegrals(K);

k  = 2;                     %   Method order.
nk = (k+1)*(k+2)*(k+3)/6;   %   Dimension of polynomial space.
                            %   Local nomber of dofs.
NK = nN + nE*(k-1) + nF*k*(k-1)/2 + k*(k^2-1)/6;
% 
% X = (X-repmat(Kc,size(X,1),1))/hK;

B = zeros(nk, NK);
intPos = G.cells.BintPos(K):G.cells.BintPos(K+1)-1;
dofVec = [nodes', edges' + G.nodes.num, ...
          faces' + G.nodes.num + G.edges.num];
B(2:nk,1:NK-1) = G.cells.Bint(intPos, dofVec);
B(1,NK) = 1;
B([5,8,10],NK) = -2*vol/hK.^2;

D = [monomialNodeVals; monomialEdgeVals;  ...
     monomialFaceInt ; monomialCellInt/vol]; 
 
%             m3D =      @(X) [ones(size(X,1),1) , ...
%                 (X(:,1)-Kc(1))/hK              , ...   %   (1,0,0)
%                (X(:,2)-Kc(2))/hK               , ...   %   (0,1,0)
%                (X(:,3)-Kc(3))/hK               , ...   %   (0,0,1)
%                (X(:,1)-Kc(1)).^2/hK^2          , ...   %   (2,0,0)
%                (X(:,1)-Kc(1)).*(X(:,2)-Kc(2))/hK^2 , ...   %   (1,1,0)
%                (X(:,1)-Kc(1)).*(X(:,3)-Kc(3))/hK^2, ...   %   (1,0,1)
%                (X(:,2)-Kc(2)).^2/hK^2, ...   %   (0,2,0) 
%                (X(:,2)-Kc(2)).*(X(:,3)-Kc(3))/hK^2, ...   %   (0,1,1)
%                (X(:,3)-Kc(3)).^2/hK^2];      %   (0,0,2)icenter of K.
%            
% Xm = (X-repmat(Kc,size(X,1),1))./hK;
% Ecm = (Ec-repmat(Kc,size(Ec,1),1))./hK;
% mFI= bsxfun(@rdivide, polygonInt3D(G,faces,m3D),faceAreas);
% mCI = polyhedronInt(G,K,m3D)/vol;
% D-[m3D([X;Ec]);mFI;mCI]


 
 
M = B*D;

PNstar = M\B;
PN = D*PNstar;

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

g = @(X) X(:,1).^2 + X(:,3).*X(:,2)/27*60 + 10;
g = @(X) ones(size(X,1),1);

gF = polygonInt3D(G,faces,g);
gK = polyhedronInt(G,K,g);
gv = [g([X;Ec]); gF./faceAreas; gK/vol];

er1 = max(abs((gv - PN*gv)./gv));
er2 = max(abs(g(X) - m3D(X)*PNstar*gv)./g(X));

Mtilde = [zeros(1,nk); M(2:nk,:)];

Sl = PNstar'*Mtilde*PNstar + hK*(eye(NK)-PN)'*(eye(NK)-PN);

%%  LOCAL LOAD TERM.                                                     %%

                            %   Matrix of integrals over K of all
                            %   combinations of linear monomials. 

H = [monomialCellInt([1,2,3,4])  ; ...
     monomialCellInt([2,5,6,7])  ; ...
     monomialCellInt([3,6,8,9])  ; ...
     monomialCellInt([4,7,9,10])];             
                            %   \Pi^\Nabla in the monomial basis
                            %   \mathcal{M}_1.
PNstar = M(1:4,1:4)\B(1:4,:);
                            %   Local load term.
                            
fChi = [f([X; Ec]); ...
        fFaceIntegrals./faceAreas   ; ...
        fCellIntegral/vol];

bl = PNstar'*H*PNstar*fChi;

dofVec = [nodes', edges' + G.nodes.num, ...
          faces' + G.nodes.num + G.edges.num, ...
          K + G.nodes.num + G.edges.num + G.faces.num];
      
end