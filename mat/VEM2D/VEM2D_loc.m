function [AP, bP, dofVec, PNstar] = VEM2D_loc(G, P, KP, mu, rho, k, ...
                               rate, f, sigma, cartGridQ, m, grad_m, int_m)
%   Calculates local stiffness matrix and load term for the virtual element
%   method for the 2D Poisson equation.
%
%   SYNOPSIS:
%       [AK, bK, dofVec, PNstar] = VEM2D_loc(G, f, K, k, alpha, rate)
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
%       [1]     - Ø. S. Klemetsdal: 'The virtual element method as a common
%                 framework for finite element and finite difference
%                 methods - Numerical and theoretical analysis'. MA thesis.
%                 Norwegian University of Science and Technology.
%-----------------------------------------------------------------ØSK-2016-

%{
   Copyright (C) 2016 Øystein Strengehagen Klemetsdal. See Copyright.txt
   for details.
%}

%%  CELL DATA                                                            %%

Pc = G.cells.centroids(P,:);
hP = G.cells.diameters(P);
aP = G.cells.volumes(P);
KP = [KP([1,3]);KP([2,4])]; % Permeability tensor, transposed.

edgeNum = G.cells.facePos(P):G.cells.facePos(P+1)-1;
edges   = G.cells.faces(edgeNum);
if size(edges,1) == 1;
    edges = edges';
end
nE          = size(edges,1);
Ec          = G.faces.centroids(edges,:);
hE          = G.faces.areas(edges,:);
edgeNormals = G.faces.normals(edges,:);
edgeSign    = (-ones(nE,1)).^(G.faces.neighbors(edges,1) ~= P); 
edgeNormals = bsxfun(@times, edgeNormals, edgeSign);

nodeNum = mcolon(G.faces.nodePos(edges),G.faces.nodePos(edges+1)-1);
nodes   = G.faces.nodes(nodeNum);
if size(nodes,1) == 1
    nodes = nodes';
end
nodes   = reshape(nodes,2,[])';
nN      = size(nodes,1);
nodes(edgeSign == -1,:) = nodes(edgeSign == -1,2:-1:1);
nodes   = nodes(:,1);

nk = (k+1)*(k+2)/2;
nkk = k*(k+1)/2;
NK = nN + nE*(k-1) + k*(k-1)/2;

X = [G.nodes.coords(nodes,:); Ec];

Xmon = bsxfun(@minus,X,Pc)/hP;                     

%%  CALCULATE B AND D MATRICES                                           %%

if k == 1
    
    D = m(Xmon(1:nN,:));
    
    XintB = repmat(Xmon(nN+1:end,:),2,1);
    edgeNormals = [edgeNormals;...
                   [edgeNormals(nE,:);edgeNormals(1:nE-1,:)]];
    
    intB = .5*sum(grad_m(XintB)*KP.*repmat(edgeNormals,2,1),2);
    intB = reshape(intB,2*nE,2);
    intB = (intB(1:nE,:) + intB(nE+1:2*nE,:))/hP;
    
%     B = [ones(1,NK)/NK; intB'];
    B = [.5*(hE(1:end)+hE([end,1:end-1]))'; intB'];
    
    H = aP;
    if isa(f,'function_handle')
        fHat = f(X(1:nN,:));
    else
        fHat = f*ones(nN,1);
    end
    fHat = mu*(fHat + rate/aP);
    rateVec = 0;
    
    dofVec = nodes;
    
elseif k == 2
    
    intD = bsxfun(@times                               , ...
                  (int_m(Xmon(1:nN,:))                          ...
                   + [int_m(Xmon(2:nN,:));int_m(Xmon(1,:))])/6  ...
                   +  int_m(Xmon(nN+1:end,:))*2/3              , ...
                   edgeNormals(:,1));
    intD = sum(intD,1)*hP/aP;
    D = [m(Xmon); intD];

    XintB = [Xmon(1:nN,:); Xmon(1:nN,:); Xmon(nN+1:end,:)];
    edgeNormals = [edgeNormals;...
                   [edgeNormals(nE,:);edgeNormals(1:nE-1,:)]; ...
                   edgeNormals];
    
    intB = sum(grad_m(XintB)*KP.*repmat(edgeNormals,5,1),2);
    intB = reshape(intB,3*nN,5);
    intB = [(intB(1:nN,:) + intB(nN+1:2*nN,:))/6; ...
             intB(2*nN+1:end,:)*2/3]/hP;
         
    B = zeros(nk, NK);
    B(1,NK) = 1;
    B(2:nk, 1:NK-1) = intB';
    B(4,NK) = B(4,NK) - 2*aP*KP(1,1)/hP^2;
    B(6,NK) = B(6,NK) - 2*aP*KP(2,2)/hP^2;
%     B([4, 6], NK) = B([4, 6], NK) - 2*aP/hP^2;

    H = zeros(nkk, nkk);
    H(1,:) = intD(:, [1,2,3])*aP;
    H(2,:) = intD(:, [2,4,5])*aP;
    H(3,:) = intD(:, [3,5,6])*aP;
    
    if isa(f,'function_handle')
        fInt = polygonInt(G, P, f, k+1)/aP;
        fHat = mu*[f(X); fInt];
    else
        fHat = mu*f*ones(2*nN+1,1);
    end
    rateVec = zeros(NK,1);
    rateVec(NK) = mu*rate;
    
    dofVec = [nodes', edges' + G.nodes.num, P + G.nodes.num + G.faces.num];
    
end

%%  CALCULATE LOCAL STIFFNESS MATRIX AND LOAD TERM                       %%


M = B*D;

PNstar = M\B;
PN = D*PNstar;
Mtilde = [zeros(1,nk) ; M(2:nk,:)];

if cartGridQ
    hx = abs(max(X(:,1))-min(X(:,1)));
    hy = abs(max(X(:,2))-min(X(:,2)));
    Q = sqrt(9/(4*hx*hy))*[-1,1,-1,1]';
else
    Q = orth(eye(NK) - PN);
end

sigma = diag(sigma,0);
AP = PNstar'*Mtilde*PNstar + trace(KP)*(eye(NK)-PN)'*Q*sigma*Q'*(eye(NK)-PN);

PNstar0 = M(1:nkk,1:nkk)\B(1:nkk,:);
bP = PNstar0'*H*PNstar0*fHat + rateVec;

end