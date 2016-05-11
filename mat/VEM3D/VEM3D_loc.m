function [AK, bK, dofVec, PNstar] ...
                           = VEM3D_loc(G, K, f, k, m, sigma, rate, mu, rho)
%--------------------------------------------------------------------------
%   Calculates local stiffness matrix and load term for the virtual element
%   method for the 2D Poisson equation.
%
%   SYNOPSIS:
%       [AK, bK, dofVec, PNstar] ...
%                             = VEM2D_loc(G, f, K, k, alpha, rate, mu, rho)
%
%   DESCRIPTION:
%       Calculates local stiffness matrix and load term for cell K of grid
%       G for the virtual element method of order k for the Poisson
%       equation
%
%           -\Delta u = f,
%
%       or, if a fluid is specified,
%
%           -\Delta p = \frac{\mu}{\rho} f,
%
%   REQUIRED PARAMETERS:
%       G       - MRST grid.
%       K       - Cell to build local stiffness matrix for.
%       f       - Source term. Either a function handle, or a scalar. In
%                 the latter case it is interpreted as a constant function.
%       m       - Monomials, see funtion retrieveMonomials.
%       grad_m  - Monomial gradients, see funtion retrieveMonomials.
%       int_m   - Monomial anti-derivatives, see funtion retrieveMonomials.
%       k       - Method order. Supported orders are k = 1 and k = 2.
%       sigma   - Vector of nker constant for scaling of the local load
%                 term. See [1] for detials.
%       rate    - Production rate in K. Set to 0 if K is not a source cell.
%       mu      - Dynamic viscosity of fluid. Set to 1 if N/A
%       rho     - Fluid density. Set to 1 if N/A.
%
%   RETURNS:
%       AK      - Local stiffness matrix.
%       bK      - Local load term.
%       dofVec  - Map from local to global stiffness matrix.
%       PNstar  - Projection operator \Pi^\nabla in the monomial basis
%                 \mathcal{M]_k(K). See [1] for details.
%
%   REFERENCES:
%       [1]     - Thesis title.
%-----------------------------------------------------------------ØSK-2016-

%{
   Copyright (C) 2016 Øystein Strengehagen Klemetsdal. See Copyright.txt
   for details.
%}

%%  CELL DATA                                                            %%

                            %   Node data for cell K.
nodeNum = G.cells.nodePos(K):G.cells.nodePos(K+1)-1;
nodes   = G.cells.nodes(nodeNum);
if size(nodes,1) == 1;
    nodes = nodes';
end
if k == 1
    X = G.nodes.coords(nodes,:);
end
nN      = size(nodes,1);
Kc  = G.cells.centroids(K,:);
hK  = G.cells.diameters(K);
aK = G.cells.volumes(K);
nE = 0; nF = 0;

faceNum     = G.cells.facePos(K):G.cells.facePos(K+1)-1;
faces       = G.cells.faces(faceNum);
if size(faces,1) == 1;
    faces = faces';
end

if k == 2
                                %   Edge data for cell K.
    edgeNum = G.cells.edgePos(K):G.cells.edgePos(K+1)-1;
    edges   = G.cells.edges(edgeNum);
    if size(edges,1) == 1;
        edges = edges';
    end
    nE      = size(edges,1);

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
    fCellIntegral = G.cells.fCellIntegrals(K);
    
    X = [G.nodes.coords(nodes,:); G.edges.centroids(edges,:)];
        
end

Xmon = bsxfun(@minus, X, Kc)/hK;

nk  = (k+1)*(k+2)*(k+3)/6;
nkk = k*(k+1)*(k+2)/6;
NK  = nN + nE*(k-1) + nF*k*(k-1)/2 + k*(k^2-1)/6;

%%  BUILD MATRIX B AND D                                                 %%

B = zeros(nk, NK);
intPos = G.cells.BintPos(K):G.cells.BintPos(K+1)-1;

% m3D = @(X) [ones(size(X,1),1), ...
%             (X(:,1)-Kc(1))/hK, ...
%             (X(:,2)-Kc(2))/hK, ...
%             (X(:,3)-Kc(3))/hK, ...
%             (X(:,1)-Kc(1)).^2/hK^2, ...
%             (X(:,1)-Kc(1)).*(X(:,2)-Kc(2))/hK^2, ...
%             (X(:,1)-Kc(1)).*(X(:,3)-Kc(3))/hK^2, ...
%             (X(:,2)-Kc(2)).^2/hK^2, ...
%             (X(:,2)-Kc(2)).*(X(:,3)-Kc(3))/hK^2, ...
%             (X(:,3)-Kc(3)).^2/hK^2];
    
% cellIntegrals = polyhedronInt(G, K, m3D, 2);

if k == 1

%     B(1,:) = 1/nFF;      % CHECK!
    
    dofVec = nodes';
    B(1,:) = sum(G.faces.B1int(faces,dofVec),1)/sum(G.faces.areas(faces));
    B(2:nk,:) = G.cells.Bint(intPos, dofVec);
    
    D = m(Xmon);

    H = aK;
    if isa(f,'function_handle')
        fHat = mu/rho*f(X);
    else
        fHat = mu/rho*f*ones(nN,1);
    end
    fHat = fHat + mu/rho*rate/aK;
    rateVec = 0;
    
    dofVec = nodes;
    
elseif k == 2
    B(1,NK) = 1;
    B(2:4, nN + nE*(k-1) + 1: nN + nE*(k-1) + nF*k*(k-1)/2) = ...
    faceNormals'/hK;
    dofVec = [nodes', edges' + G.nodes.num, ...
              faces' + G.nodes.num + G.edges.num];
    B(5:nk,1:NK-1) = G.cells.Bint(intPos, dofVec);
    B([5,8,10],NK) = -2*aK/hK.^2;
    
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
    
    faceIntegrals = polygonInt3D(G, faces, m3D, 2);
    cellIntegrals = polyhedronInt(G, K, m3D, 2);
    
    D = [m(Xmon)                                    ; ...
         bsxfun(@rdivide, faceIntegrals, faceAreas) ; ...
         cellIntegrals/aK                          ];

    H = [cellIntegrals([1,2,3,4])  ; ...
         cellIntegrals([2,5,6,7])  ; ...
         cellIntegrals([3,6,8,9])  ; ...
         cellIntegrals([4,7,9,10])];             

    if isa(f,'function_handle')
        fHat = mu/rho*[f(X); ...
                       fFaceIntegrals./faceAreas;
                       fCellIntegral/aK];
    else
        fHat = mu/rho*f*ones(nN + nE + nF +1,1);
    end
    rateVec = zeros(NK,1);
    rateVec(NK) = mu/rho*rate;
    
    dofVec = [nodes', edges' + G.nodes.num, ...
              faces' + G.nodes.num + G.edges.num, ...
              K + G.nodes.num + G.edges.num + G.faces.num];

end

%%  LOCAL STIFFNESS MATRIX                                               %%
 
M = B*D;
PNstar = M\B;
PN = D*PNstar;

Mtilde = [zeros(1,nk); M(2:nk,:)];

Q = orth(eye(NK)-PN);

sigma = diag(sigma,0);
AK = PNstar'*Mtilde*PNstar + hK*(eye(NK)-PN)'*Q*sigma*Q'*(eye(NK)-PN);

PNstar0 = M(1:nkk,1:nkk)\B(1:nkk,:);

bK = PNstar0'*H*PNstar0*fHat + rateVec;

% PNstar0 = M(1:nkk,1:nkk)\B(1:nkk,:);
% if k == 1
%     Pf = polyhedronInt(G,K,f,7)/aK;
%     bK = PNstar0'*H*Pf + rateVec;
% elseif k == 2
%     bK = PNstar0'*H*PNstar0*fHat + rateVec;
% end
% elseif k == 1
%     
%     [Xq, w, V, vol] = tetrahedronQuadRule(7);
%     
%     Vdiff = V(1:end-1,:) - V(2:end,:);
%     nq = size(Xq,1);
%     
%     m = retrieveMonomials(3,k);
%     tri = delaunay(X);
%     nTri = size(tri,1);
%     
%     X = X(reshape(tri',[],1),:);
%     ii = mcolon(1:4:4*nTri, 3:4:4*nTri);
%     Xdiff = X(ii,:) - X(ii+1,:);
% 
%     ii = repmat((1:3)',3*nTri,1);
%     add = repmat((0:1:nTri-1)*3,9,1);
%     ii = ii + add(:);
%     jj = repmat((1:3), 3, 1);
%     jj = repmat(jj(:), nTri,1);
%     jj = jj(:) + add(:);
%     Vd = sparse(ii,jj,repmat(Vdiff(:),nTri,1), 3*nTri, 3*nTri);
%     Xdiff = Xdiff';
%     Xd = sparse(jj,ii,Xdiff(:), 3*nTri, 3*nTri);
%     
%     A = Vd\Xd;
% 
%     Xb = reshape(X(1:4:end,:)',1,[]);
% 
%     b = Xb - repmat(V(1,:),1,nTri)*A;
% 
%     [n,~] = size(A);
%     
%     main = ones(n,1);
%     sub1 = main;
%     sub1(3:3:end) = 0;
%     sub2 = zeros(n,1);
%     sub2(1:3:end) = 1;
%     Id = spdiags([sub2, sub1, main,sub1(end:-1:1),sub2(end:-1:1)], -2:1:2, n,n);
% 
%     Ad = full(A(Id == 1));
%     Ad = reshape(Ad,3,3,[]);
%     D = squeeze(Ad(1,1,:).*Ad(2,2,:).*Ad(3,3,:) + ...
%                 Ad(2,1,:).*Ad(3,2,:).*Ad(1,3,:) + ...
%                 Ad(3,1,:).*Ad(1,2,:).*Ad(2,3,:) - ...
%                 Ad(3,1,:).*Ad(2,2,:).*Ad(1,3,:) - ...
%                 Ad(2,1,:).*Ad(1,2,:).*Ad(3,3,:) - ...
%                 Ad(1,1,:).*Ad(3,2,:).*Ad(2,3,:));
%     D = abs(D);
%     
%     XhatTmp = bsxfun(@plus, repmat(Xq, 1, nTri)*A, b);
%     Xhat = zeros(nq*nTri,3);
%     Xhat(:,1) = reshape(XhatTmp(:,1:3:3*nTri),[],1);
%     Xhat(:,2) = reshape(XhatTmp(:,2:3:3*nTri),[],1);
%     Xhat(:,3) = reshape(XhatTmp(:,3:3:3*nTri),[],1);
%     
%     Xmon = bsxfun(@minus, Xhat,Kc)/hK;
%     PNphi = m(Xmon)*PNstar;
%     
%     bK = (vol*repmat(w,1,nTri).*rldecode(D,nq*ones(nTri,1),1)'*bsxfun(@times,PNphi,f(Xhat)))';
% end
        
% % 
% g = @(X) X(:,1).^5 + X(:,3).*X(:,2);
% PNgInt = aK*PNstar0*g(X);
% gInt = polyhedronInt(G,K,g,7);
% if abs(PNgInt-gInt) > 10e-10
%     a = abs(PNgInt-gInt)/gInt
%     sum(PNstar0)
% end

end