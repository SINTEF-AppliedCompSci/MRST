function [A, b, PNstarT] ...
                 = VEM3D_glob(G, f, k, bc, sigma, projectors, src, mu, rho)
%--------------------------------------------------------------------------
%   Assmebles the global stiffness matrix and load term for the virtual
%   element method for the 2D Poisson equation.
%
%   SYNOPSIS:
%       [A, b, PNstarT] ...
%                = VEM2D_glob(G, f, k, bc, alpha, projectors, src, mu, rho)
%
%   DESCRIPTION:
%       Assmebles the global stiffness matrix and load term for the virtual
%       element method of order k on grid G for the 2D Poisson equation 
%
%           -\Delta u = f,
%
%       or, if a fluid is specified,
%
%           -\Delta p = \frac{\mu}{\rho} f,
%
%   REQUIRED PARAMETERS:
%       G          - MRST grid.
%       f          - Source term. Either a function handle, or a scalar. In
%                    the latter case it is interpreted as a constant
%                    function.
%       k          - Method order. Supported orders are k = 1 and k = 2.
%       bc         - Struct of boundary conditions constructed using
%                    VEM2D_addBC.
%       sigma      - G.cells.num x nker matrix of constants for scaling
%                    of the local load terms.
%                    nker = \dim \ker \Pi^\nabla. See [1] for detail.
%       projectors - Boolean. If true, matrix representations
%                    of \Pi^\nabla in the monomial basis \mathcal_k(K) will
%                    be stored.
%       src        - Source term struct constructed using addSource.
%       mu         - Dynamic viscosity of fluid. Set to 1 if N/A
%       rho        - Fluid density. Set to 1 if N/A.
%
%   RETURNS:
%       A          - Global stiffness matrix.
%       b          - Global load term.
%       PNstarT    - Projection operators \Pi^\nabla in the monomial basis
%                    \mathcal{M]_k(K). Stored as transpose. See [1] for
%                    details.
%
%   REFERENCES:
%       [1]     - Thesis title.
%-----------------------------------------------------------------ØSK-2016-

%{
   Copyright (C) 2016 Øystein Strengehagen Klemetsdal. See Copyright.txt
   for details.
%}

%%  SET FUNCTION SPACE DIMS                                              %%

nN = G.nodes.num;
nE = G.edges.num;
nF = G.faces.num;
nK = G.cells.num;

nk   = (k+1)*(k+2)*(k+3)/6;
NK   = diff(G.cells.nodePos) + diff(G.cells.edgePos)*(k-1) + ...
       diff(G.cells.facePos)*k*(k-1)/2 + k*(k^2-1)/6;
nker = NK - nk;
N  = nN + nE*(k-1) + nF*k*(k-1)/2 + nK*k*(k^2-1)/6;

%%  COMPUTE FACE PROJECTORS, MONOMIAL DOFS AND SOURCE TERM INTEGRALS     %%

fprintf('... computing face projectors\n');

I = faceProjectors(G,k);
nk2D = (k+1)*(k+2)/2;
BintPos = (0:nk2D:nk2D*G.cells.num) + 1;
G.cells.('Bint') = I;
G.cells.('BintPos') = BintPos;

fprintf('... computing monomial values\n');

monomialVals = monomialValues(G,k);
monomialNodeValsPos = [1, cumsum(diff(G.cells.nodePos)') + 1];
G.cells.('monomialNodeVals') = monomialVals(1:monomialNodeValsPos(end)-1,:);
G.cells.('monomialNodeValsPos') = monomialNodeValsPos;

if k == 2
    monomialEdgeValsPos = [1, cumsum(diff(G.cells.edgePos)') + 1];
    G.cells.('monomialEdgeVals') = monomialVals(monomialNodeValsPos(end):end,:);
    G.cells.('monomialEdgeValsPos') = monomialEdgeValsPos;
    
    fprintf('... computing source term integrals\n');
    
    IFf = polygonInt3D(G,1:G.faces.num,f,k+1);
    G.faces.('fFaceIntegrals') = IFf;
    ICf = polyhedronInt(G,1:G.cells.num,f,k+1);
    G.cells.('fCellIntegrals') = ICf;
end

%%  CREATE MAPPINGS                                                      %%

dofPosA = [1, cumsum((NK.^2)') + 1];
dofPosb = [1, cumsum( NK'    ) + 1];

%%  ADJUST INPUT AND OUTPUT PARAMETRES                                   %%

if numel(sigma) > 1
    sigmaPos = [1, cumsum(nker') + 1];
end

rate = zeros(G.cells.num,1);
if ~isempty(src)
    rate(src.cell) = src.rate;
end

if projectors
    PNstarT = zeros(dofPosb(end)-1,nk);
else
    PNstarT = 0;
end

%%  ASSEMBLE GLOBAL STIFFNESS MATRIX AND LOAD VECTOR                     %%

step = floor(nK/10);

iiA = zeros(1,dofPosA(end)-1);
jjA = zeros(1,dofPosA(end)-1);
AVec = zeros(1,dofPosA(end)-1);

iib = zeros(1,dofPosb(end)-1);
bVec = zeros(1,dofPosb(end)-1);


fprintf('Computing local block matrices ...\n')
tic;

for K = 1:nK

    if rem(K,step) == 0
        fprintf('... Calculating local block matrix for cell %d\n', K);
    end 
    
    if numel(sigma) > 1
        sigmaK = sigma(sigmaPos(K):sigmaPos(K+1)-1);
    else
        sigmaK = sigma;
    end
    
    [AK, bK, dofVec, PNstar] ...
       = VEM3D_loc(G, K, f, k, sigmaK, rate(K), mu, rho);

    NK = numel(dofVec);
    
    iiK = repmat(dofVec', NK, 1);
    jjK = repmat(dofVec , NK, 1);
    jjK = jjK(:);
    
    iiA(dofPosA(K):dofPosA(K+1)-1) = iiK;
    jjA(dofPosA(K):dofPosA(K+1)-1) = jjK;
    AVec(dofPosA(K):dofPosA(K+1)-1)= AK(:);

    iib(dofPosb(K):dofPosb(K+1)-1) = dofVec;
    bVec(dofPosb(K):dofPosb(K+1)-1) = bK;
    
    if projectors
        PNstarT(dofPosb(K):dofPosb(K+1)-1,:) = PNstar';
    end
    
end

stop = toc;
fprintf('Done in %f seconds.\n\n', stop);

A = sparse(iiA, jjA, AVec, N, N);
b = sparse(iib, ones(1, numel(iib)), bVec, N, 1);

%%  APPLY BOUNDARY CONDITIONS                                            %%

[A,b] = VEM3D_bc(G,A,b,bc,k);

end


% %% OLD
% nN = G.nodes.num;
% nE = G.edges.num;
% nF = G.faces.num;
% nK = G.cells.num;
% 
% N = nN + nE*(k-1) + nF*k*(k-1)/2 + nK*k*(k^2-1)/6;
% 
% if k == 1
%     
%     dofPosA = [1, cumsum((diff(G.cells.nodePos')).^2) + 1];
%     dofPosb = [1, cumsum(diff(G.cells.nodePos')) + 1];
%         
% elseif k == 2
% 
%     dofPosA = [1, cumsum((diff(G.cells.nodePos') + ...
%                           diff(G.cells.edgePos') + ...
%                           diff(G.cells.facePos') + 1).^2) + 1];
%     dofPosb = [1, cumsum(diff(G.cells.nodePos') + ...
%                          diff(G.cells.edgePos') + ...
%                          diff(G.cells.facePos') + 1) + 1];
% 
% end
%                      
% step = floor(nK/10);
% 
% % itS = 1; itb = 1;
% 
% iiA = zeros(1,dofPosA(end)-1);
% jjA = zeros(1,dofPosA(end)-1);
% 
% AVec = zeros(1,dofPosA(end)-1);
% bVec = zeros(1,dofPosb(end)-1);
% 
% % fprintf('Computing local block matrices ...\n')
% % tic;
% 
% for K = 1:nK
% 
%     if rem(K,step) == 0
%         fprintf('... Calculating local block matrix for cell %d\n', K);
%     end
%     
%     [AK, bK, dofVec] = VEM3D_loc(G, f, K, k);
%     
%     if K == 1
%         G.cells.('AK') = AK;
%     end
% 
%     NK = numel(dofVec);
%     
%     iiK = repmat(dofVec', NK, 1);
%     jjK = repmat(dofVec , NK, 1);
%     jjK = jjK(:);
%     
%     iiA(dofPosA(K):dofPosA(K+1)-1) = iiK;
%     jjA(dofPosA(K):dofPosA(K+1)-1) = jjK;
%     AVec(dofPosA(K):dofPosA(K+1)-1)= AK(:);
% 
%     iib(dofPosb(K):dofPosb(K+1)-1) = dofVec;
%     bVec(dofPosb(K):dofPosb(K+1)-1) = bK(:);
%     
% end
% 
% % stop = toc;
% % fprintf('Done in %f seconds.\n\n', stop);
% 
% S = sparse(iiA, jjA, AVec, N, N);
% b = sparse(iib, ones(1, numel(iib)), bVec);
% 
% [bcDof, bBC] = VEM3D_bc(G, bc, k);
% b(bcDof == 1) = bBC(bcDof == 1);
% SBC = spdiags(ones(N,1),0,N,N);
% S(bcDof == 1,:) = SBC(bcDof == 1,:);
% b(bcDof == 2) = b(bcDof == 2) + bBC(bcDof == 2);