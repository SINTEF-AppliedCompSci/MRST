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
nF          = size(faces,1);
faceAreas   = G.faces.areas(faces);
faceNormals = G.faces.normals(faces,:);
faceSigns    = (-ones(nF,1)).^(G.faces.neighbors(faces,1) ~= K);
faceNormals = bsxfun(@times, faceNormals,faceSigns);
fFaceIntegrals = G.faces.fFaceIntegrals(faces);

                            %   Cell data for cell K.
Kc  = G.cells.centroids(K,:);
hK  = G.cells.diameters(K);
vol = G.cells.volumes(K);
fCellIntegral = G.cells.fCellIntegrals(K);


k  = 2;                     %   Method order.
nk = (k+1)*(k+2)*(k+3)/6;   %   Dimension of polynomial space.
                            %   Local nomber of dofs.
NK = nN + nE*(k-1) + nF*k*(k-1)/2 + k*(k^2-1)/6;


%%  BUILD MATRIX B                                                       %%

B = zeros(nk, NK);
intPos = G.cells.BintPos(K):G.cells.BintPos(K+1)-1;
dofVec = [nodes', edges' + G.nodes.num, ...
          faces' + G.nodes.num + G.edges.num];
% full(G.cells.Bint(intPos, :));
B(5:nk,1:NK-1) = G.cells.Bint(intPos, dofVec);
B(1,NK) = 1;
B(2:4, nN + nE*(k-1) + 1: nN + nE*(k-1) + nF*k*(k-1)/2) = ...
    faceNormals'/hK;
B([5,8,10],NK) = -2*vol/hK.^2;

%%  BUILD MATRIX D                                                       %%

m3D = @(X) [ones(size(X,1),1), ...
                (X(:,1)-Kc(1))/hK, ...
                (X(:,2)-Kc(2))/hK, ...
                (X(:,3)-Kc(3))/hK, ...
                (X(:,1)-Kc(1)).^2/hK^2, ...
                (X(:,1)-Kc(1)).*(X(:,2)-Kc(2))/hK^2, ...
                (X(:,1)-Kc(1)).*(X(:,3)-Kc(3))/hK^2, ...
                (X(:,2)-Kc(2)).^2/hK^2, ...
                (X(:,2)-Kc(2)).*(X(:,3)-Kc(3))/hK^2, ...
                (X(:,3)-Kc(3)).^2/hK^2];

faceIntegrals = polygonInt3D(G, faces, m3D);
cellIntegrals = polyhedronInt(G, K, m3D);
            
D = [monomialNodeVals                           ; ...
     monomialEdgeVals                           ; ...
     bsxfun(@rdivide, faceIntegrals, faceAreas) ; ...
     cellIntegrals/vol                          ];

%%  LOCAL STIFFNESS MATRIX                                               %%
 
M = B*D;
PNstar = M\B;
PN = D*PNstar;

Mtilde = [zeros(1,nk); M(2:nk,:)];

Sl = PNstar'*Mtilde*PNstar + hK*(eye(NK)-PN)'*(eye(NK)-PN);

% SK = hK*(eye(NK)-PN)'*(eye(NK)-PN);
% 
% norm(SK,'fro')

%%  LOCAL LOAD TERM.                                                     %%

                            %   Matrix of integrals over K of all
                            %   combinations of linear monomials. 
H = [cellIntegrals([1,2,3,4])  ; ...
     cellIntegrals([2,5,6,7])  ; ...
     cellIntegrals([3,6,8,9])  ; ...
     cellIntegrals([4,7,9,10])];             
                            %   \Pi^\Nabla in the monomial basis
                            %   \mathcal{M}_1.
PNstar = M(1:4,1:4)\B(1:4,:);
                            %   Dofs \chi^i(f).
fChi = [f([X; Ec]); ...
        fFaceIntegrals./faceAreas   ; ...
        fCellIntegral/vol];

bl = PNstar'*H*PNstar*fChi;


%%  LOCAL TO GLOBAL MAP                                                  %%

dofVec = [nodes', edges' + G.nodes.num, ...
          faces' + G.nodes.num + G.edges.num, ...
          K + G.nodes.num + G.edges.num + G.faces.num]; 

% Mdb = zeros(nk,nk);
% Mdb(1,:) = cellIntegrals/vol;
% %    1 2 3 4 5   6  7  8   9  10
% %   {1,x,y,z,x^2,xy,xz,y^2,yz,z^2}
% mm = zeros(nk,nk);
% mm(2,2) = vol; mm(3,3) = vol; mm(4,4) = vol;
% mm(2,[5,6,7])  = cellIntegrals([2,3,4]).*[2,1,1];
% mm(3,[6,8,9])  = cellIntegrals([2,3,4]).*[1,2,1];
% mm(4,[7,9,10]) = cellIntegrals([2,3,4]).*[1,1,2];
% mm(5,[5,6,7])  = cellIntegrals([5,6,7]).*[4,2,2]; 
% mm(6,6)        = cellIntegrals(5) + cellIntegrals(8);
% mm(6,[7,8,9])  = cellIntegrals([9,6,7]).*[1,2,1];
% mm(7,7)        = cellIntegrals(5) + cellIntegrals(10);
% mm(7,[9,10])   = cellIntegrals([6,7]).*[1,2];
% mm(8,[8,9])    = cellIntegrals([8,9]).*[4,2];
% mm(9,9)        = cellIntegrals(8) + cellIntegrals(10);
% mm(9,10)       = cellIntegrals(9)*2;
% mm(10,10)      = cellIntegrals(10)*4;

% 
% mm = mm/hK^2;
% 
% Mdb(2:nk,2:nk) = triu(mm(2:nk,2:nk)) + tril(mm(2:nk,2:nk)',-1);
% 
% norm(Mdb-M,'fro')
% (abs(Mdb-M)>10e-6)
% 
% 
% Fc = G.faces.centroids(faces,:);
% edgeNum = mcolon(G.faces.edgePos(faces),G.faces.edgePos(faces+1)-1);
% edges = G.faces.edges(edgeNum);
% edgeNormals = G.faces.edgeNormals(edgeNum,:);
% clf;
% plotGrid(G,K)
% hold on
% plot3(X(:,1), X(:,2), X(:,3),'*')
% plot3(Ec(:,1), Ec(:,2), Ec(:,3), 'o')
% Ec = G.edges.centroids(edges,:);
% % plot3(Ec(:,1) + edgeNormals(:,1), ...
% %       Ec(:,2) + edgeNormals(:,2), ...
% %       Ec(:,3) + edgeNormals(:,3), '+')
% 
% plot3(Fc(:,1),Fc(:,2), Fc(:,3), 'o')
% plot3(Fc(:,1) + faceNormals(:,1), Fc(:,2) ...
%               + faceNormals(:,2), Fc(:,3) ...
%               + faceNormals(:,3), '+')
% 
% hold off
% axis equal
      
end
     
%%  DEBUG

% g = @(X) X(:,1).^2 + X(:,3).*X(:,2)*1000/3 + 10;
% g = @(X) ones(size(X,1),1);
% 
% gF = polygonInt3D(G,faces,g);
% gK = polyhedronInt(G,K,g);
% gv = [g([X;Ec]); gF./faceAreas; gK/vol];
% 
% er1 = max(abs((gv - PN*gv)./gv))
% er2 = max(abs(g(X) - m3D(X)*PNstar*gv)./g(X))


% D = [monomialNodeVals; monomialEdgeVals;  ...
%      monomialFaceInt ; monomialCellInt/vol]; 
%  
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


% m3D =      @(X) [ones(size(X,1),1) , ...
%                 X(:,1)              , ...   %   (1,0,0)
%                X(:,2)               , ...   %   (0,1,0)
%                X(:,3)               , ...   %   (0,0,1)
%                X(:,1).^2          , ...   %   (2,0,0)
%                X(:,1).*X(:,2) , ...   %   (1,1,0)
%                X(:,1).*X(:,3), ...   %   (1,0,1)
%                X(:,2).^2, ...   %   (0,2,0) 
%                X(:,2).*X(:,3), ...   %   (0,1,1)
%                X(:,3).^2];      %   (0,0,2)icenter of K.
% 
% 
% Mdb = zeros(nk,nk);
% Mdb(1,:) = cellIntegrals/vol;
% %    1 2 3 4 5   6  7  8   9  10
% %   {1,x,y,z,x^2,xy,xz,y^2,yz,z^2}
% mm = zeros(nk,nk);
% mm(2,2) = vol/(hK^2); mm(3,3) = vol/(hK^2); mm(4,4) = vol/(hK^2);
% mm(2,5) = cellIntegrals(2)*2/(hK^2);
% mm(2,[6,7]) = cellIntegrals([3,4])/hK^2;
% mm(3,[6,8,9]) = cellIntegrals([2,3,4]).*[1/hK^2,2/hK^2, 1/hK^2];
% mm(4,[7,9,10]) = cellIntegrals([2,3,4]).*[1/hK^2,2/hK^2, 1/hK^2];
% mm(5,[5,6,7]) = cellIntegrals([5,6,7]).*[4/hK^2,2/hK^2, 2/hK^2]; 
% mm(6,6) = (cellIntegrals(5) + cellIntegrals(8))/hK^2;
% mm(6,[7,8,9]) = cellIntegrals([9,6,7])/hK^2;
% mm(7,7) = (cellIntegrals(5) + cellIntegrals(10))/hK^2;
% mm(7,[9,10]) = cellIntegrals([6,7]).*[1/hK^2,2/hK^2];
% mm(8,[8,9]) = cellIntegrals([8,9]).*[4/hK^2,2/hK^2];
% mm(9,9) = (cellIntegrals(8) + cellIntegrals(10))/hK^2;
% mm(9,10) = cellIntegrals(9)*2/hK^2;
% mm(10,10) = cellIntegrals(10)*4/hK^2;
% 
% Mdb(2:nk,2:nk) = triu(mm(2:nk,2:nk)) + tril(mm(2:nk,2:nk)',-1);
% 
% norm(Mdb-M,'fro')
% Mdb-M

%%  DELETED STUFF

% faceIntNum = G.cells.faceIntPos(K):G.cells.faceIntPos(K+1)-1;
% monomialFaceInt     = bsxfun(@rdivide, ...
%                      G.cells.monomialFaceIntegrals(faceIntNum,:), ...
%                      faceAreas);
% mIint = G.cells.monomialFaceIntegrals(faceIntNum,:)


% monomialCellInt = G.cells.monomialCellIntegrals(K,:);