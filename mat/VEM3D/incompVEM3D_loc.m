function [AP, bP, dofVec, PNstar] ...
      = incompVEM3D_loc(G, P, KP, mu, rho, k, rate, f, sigma, cartGridQ, m)
%--------------------------------------------------------------------------
%   Calculates local stiffness matrix and load term for the virtual element
%   method for the 2D Poisson equation.
%
%   SYNOPSIS:
%       [AP, bP, dofVec, PNstar] = VEM2D_loc(G, f, K, k, alpha, rate)
%
%   DESCRIPTION:
%       Calculates local stiffness matrix and load term for cell K of grid
%       G for the virtual element method of order k for the Poisson
%       equation
%
%           -\Delta u = f.
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
%
%   RETURNS:
%       AK      - Local stiffness matrix.
%       bK      - Local load term.
%       dofVec  - Map from local to global stiffness matrix.
%       PNstar  - Projection operator \Pi^\nabla in the monomial basis
%                 \mathcal{M]_k(K). See [1] for details.
%
%   REFERENCES:
%       [1] - Ø. S. Klemetsdal: 'The virtual element method as a common
%             framework for finite element and finite difference methods -
%             Numerical and theoretical analysis'. MA thesis. Norwegian
%             University of Science and Technology.
%-----------------------------------------------------------------ØSK-2016-

%{
   Copyright (C) 2016 Øystein Strengehagen Klemetsdal. See COPYRIGHT.txt
   for details.
%}

%%  CELL DATA                                                            %%

                            %   Node data for cell P.
nodeNum = G.cells.nodePos(P):G.cells.nodePos(P+1)-1;
nodes   = G.cells.nodes(nodeNum);
if size(nodes,1) == 1;
    nodes = nodes';
end
if k == 1
    X = G.nodes.coords(nodes,:);
end
nN      = size(nodes,1);
Pc  = G.cells.centroids(P,:);
hP  = G.cells.diameters(P);
aP = G.cells.volumes(P);
nE = 0; nF = 0;

KP = [KP(1:3); KP(4:6); KP(7:9)];

faceNum     = G.cells.facePos(P):G.cells.facePos(P+1)-1;
faces       = G.cells.faces(faceNum);
if size(faces,1) == 1;
    faces = faces';
end

if k == 2
                                %   Edge data for cell P.
    edgeNum = G.cells.edgePos(P):G.cells.edgePos(P+1)-1;
    edges   = G.cells.edges(edgeNum);
    if size(edges,1) == 1;
        edges = edges';
    end
    nE      = size(edges,1);

    faceNum     = G.cells.facePos(P):G.cells.facePos(P+1)-1;
    faces       = G.cells.faces(faceNum);
    if size(faces,1) == 1;
        faces = faces';
    end
    nF          = size(faces,1);
    faceAreas   = G.faces.areas(faces);
    faceNormals = G.faces.normals(faces,:);
    faceSigns    = (-ones(nF,1)).^(G.faces.neighbors(faces,1) ~= P);
    faceNormals = bsxfun(@times, faceNormals,faceSigns);
    fFaceIntegrals = G.faces.fFaceIntegrals(faces);

                                %   Cell data for cell P.
    fCellIntegral = G.cells.fCellIntegrals(P);
    
    X = [G.nodes.coords(nodes,:); G.edges.centroids(edges,:)];
        
end

Xmon = bsxfun(@minus, X, Pc)/hP;

nk  = (k+1)*(k+2)*(k+3)/6;
nkk = k*(k+1)*(k+2)/6;
NP  = nN + nE*(k-1) + nF*k*(k-1)/2 + k*(k^2-1)/6;

%%  BUILD MATRIX B AND D                                                 %%

B = zeros(nk, NP);
intPos = G.cells.BintPos(P):G.cells.BintPos(P+1)-1;

if k == 1
    
    dofVec = nodes';
    B(1,:) = sum(G.faces.B1int(faces,dofVec),1)/sum(G.faces.areas(faces));
    B(2:nk,:) = G.cells.Bint(intPos, dofVec);
    
    D = m(Xmon);

    H = aP;
    if isa(f,'function_handle')
        fHat = f(X);
    else
        fHat = f*ones(nN,1);
    end
    fHat = fHat + rate/aP;
    rateVec = 0;
    
    dofVec = nodes;
    
elseif k == 2
    B(1,NP) = 1;
    B(2:4, nN + nE*(k-1) + 1: nN + nE*(k-1) + nF*k*(k-1)/2) = ...
    faceNormals'/hP;
    dofVec = [nodes', edges' + G.nodes.num, ...
              faces' + G.nodes.num + G.edges.num];
    B(5:nk,1:NP-1) = G.cells.Bint(intPos, dofVec);
    B([5,8,10],NP) = -2*aP/hP.^2;
    
    m3D = @(X) [ones(size(X,1),1), ...
                (X(:,1)-Pc(1))/hP, ...
                (X(:,2)-Pc(2))/hP, ...
                (X(:,3)-Pc(3))/hP, ...
                (X(:,1)-Pc(1)).^2/hP^2, ...
                (X(:,1)-Pc(1)).*(X(:,2)-Pc(2))/hP^2, ...
                (X(:,1)-Pc(1)).*(X(:,3)-Pc(3))/hP^2, ...
                (X(:,2)-Pc(2)).^2/hP^2, ...
                (X(:,2)-Pc(2)).*(X(:,3)-Pc(3))/hP^2, ...
                (X(:,3)-Pc(3)).^2/hP^2];
    
    faceIntegrals = polygonInt3D(G, faces, m3D, 2);
    cellIntegrals = polyhedronInt(G, P, m3D, 2);
    
    D = [m(Xmon)                                    ; ...
         bsxfun(@rdivide, faceIntegrals, faceAreas) ; ...
         cellIntegrals/aP                          ];

    H = [cellIntegrals([1,2,3,4])  ; ...
         cellIntegrals([2,5,6,7])  ; ...
         cellIntegrals([3,6,8,9])  ; ...
         cellIntegrals([4,7,9,10])];             

    if isa(f,'function_handle')
        fHat = [f(X); ...
               fFaceIntegrals./faceAreas;
               fCellIntegral/aP];
    else
        fHat = f*ones(nN + nE + nF +1,1);
    end
    rateVec = zeros(NP,1);
    rateVec(NP) = rate;
    
    dofVec = [nodes', edges' + G.nodes.num, ...
              faces' + G.nodes.num + G.edges.num, ...
              P + G.nodes.num + G.edges.num + G.faces.num];

end

%%  LOCAL STIFFNESS MATRIX                                               %%
 
M = B*D;
PNstar = M\B;
PN = D*PNstar;

Mtilde = [zeros(1,nk); M(2:nk,:)];

Q = orth(eye(NP)-PN);

sigma = diag(sigma,0);
AP = PNstar'*Mtilde*PNstar ...
  + (KP(1,1) + KP(2,2) + KP(3,3))*hP*(eye(NP)-PN)'*Q*sigma*Q'*(eye(NP)-PN);

PNstar0 = M(1:nkk,1:nkk)\B(1:nkk,:);

bP = PNstar0'*H*PNstar0*fHat + rateVec;

end