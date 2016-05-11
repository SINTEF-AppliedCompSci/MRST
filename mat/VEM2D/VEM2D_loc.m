function [AK, bK, dofVec, PNstar] ...
 = VEM2D_loc(G, K, f, m, grad_m, int_m, k, sigma, cartGridQ, rate, mu, rho)
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

Kc = G.cells.centroids(K,:);
hK = G.cells.diameters(K);
aK = G.cells.volumes(K);

edgeNum = G.cells.facePos(K):G.cells.facePos(K+1)-1;
edges   = G.cells.faces(edgeNum);
if size(edges,1) == 1;
    edges = edges';
end
nE          = size(edges,1);
Ec          = G.faces.centroids(edges,:);
hE          = G.faces.areas(edges,:);
edgeNormals = G.faces.normals(edges,:);
edgeSign    = (-ones(nE,1)).^(G.faces.neighbors(edges,1) ~= K); 
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

Xmon = bsxfun(@minus,X,Kc)/hK;                     

%%  CALCULATE B AND D MATRICES                                           %%

if k == 1
    
    D = m(Xmon(1:nN,:));
    
    XintB = repmat(Xmon(nN+1:end,:),2,1);
    edgeNormals = [edgeNormals;...
                   [edgeNormals(nE,:);edgeNormals(1:nE-1,:)]];
    
    
    intB = .5*sum(grad_m(XintB).*repmat(edgeNormals,2,1),2);
    intB = reshape(intB,2*nE,2);
    intB = (intB(1:nE,:) + intB(nE+1:2*nE,:))/hK;
    B = [ones(1,NK)/NK; intB'];
    
    H = aK;
    if isa(f,'function_handle')
        fHat = mu/rho*f(X(1:nN,:));
    else
        fHat = mu/rho*f*ones(nN,1);
    end
    fHat = fHat + mu/rho*rate/aK;
    rateVec = 0;
    
    dofVec = nodes;
    
elseif k == 2
    
    intD = bsxfun(@times                               , ...
                  (int_m(Xmon(1:nN,:))                          ...
                   + [int_m(Xmon(2:nN,:));int_m(Xmon(1,:))])/6  ...
                   +  int_m(Xmon(nN+1:end,:))*2/3              , ...
                   edgeNormals(:,1));
    intD = sum(intD,1)*hK/aK;
    D = [m(Xmon); intD];

    XintB = [Xmon(1:nN,:); Xmon(1:nN,:); Xmon(nN+1:end,:)];
    edgeNormals = [edgeNormals;...
                   [edgeNormals(nE,:);edgeNormals(1:nE-1,:)]; ...
                   edgeNormals];
    
    intB = sum(grad_m(XintB).*repmat(edgeNormals,5,1),2);
    intB = reshape(intB,3*nN,5);
    intB = [(intB(1:nN,:) + intB(nN+1:2*nN,:))/6; ...
             intB(2*nN+1:end,:)*2/3]/hK;
         
    B = zeros(nk, NK);
    B(1,NK) = 1;
    B(2:nk, 1:NK-1) = intB';
    B([4, 6], NK) = B([4, 6], NK) - 2*aK/hK^2;

    H = zeros(nkk, nkk);
    H(1,:) = intD(:, [1,2,3])*aK;
    H(2,:) = intD(:, [2,4,5])*aK;
    H(3,:) = intD(:, [3,5,6])*aK;
    
    if isa(f,'function_handle')
        fInt = polygonInt(G, K, f, k+1)/aK;
        fHat = mu/rho*[f(X); fInt];
    else
        fHat = mu/rho*f*ones(2*nN+1,1);
    end
    rateVec = zeros(NK,1);
    rateVec(NK) = mu/rho*rate;
    
    dofVec = [nodes', edges' + G.nodes.num, K + G.nodes.num + G.faces.num];
    
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
AK = PNstar'*Mtilde*PNstar + (eye(NK)-PN)'*Q*sigma*Q'*(eye(NK)-PN);

PNstar0 = M(1:nkk,1:nkk)\B(1:nkk,:);
bK = PNstar0'*H*PNstar0*fHat + rateVec;

end