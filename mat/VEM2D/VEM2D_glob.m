function [A, b, PNstarT] ...
                 = VEM2D_glob(G, f, k, bc, alpha, projectors, src, mu, rho)
%--------------------------------------------------------------------------
%   Assmebles the global stiffness matrix and load term for the virtual
%   element method for the Poisson equation.
%
%   SYNOPSIS:
%       [A, b, PNstarT] ...
%                = VEM2D_glob(G, f, k, bc, alpha, projectors, src, mu, rho)
%
%   DESCRIPTION:
%       Assmebles the global stiffness matrix and load term for the virtual
%       element method of order k on grid G for the poisson equation 
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
%       alpha      - G.cells.num x 1 matrix of constants for scaling of the
%                    local load terms.
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

nN = G.nodes.num;
nE = G.faces.num;
nK = G.cells.num;

N = nN + nE*(k-1) + nK*k*(k-1)/2;

if k == 1
    dofPosA = [1, cumsum((diff(G.cells.nodePos')).^2) + 1];
    dofPosb = [1, cumsum( diff(G.cells.nodePos')) + 1];
elseif k == 2
    dofPosA = [1, cumsum((diff(G.cells.nodePos') + ...
                          diff(G.cells.facePos') + 1).^2) + 1];
    dofPosb = [1, cumsum( diff(G.cells.nodePos') + ...
                          diff(G.cells.facePos') + 1) + 1];
end
                     
step = floor(nK/10);

iiA = zeros(1,dofPosA(end)-1);
jjA = zeros(1,dofPosA(end)-1);
AVec = zeros(1,dofPosA(end)-1);

iib = zeros(1,dofPosb(end)-1);
bVec = zeros(1,dofPosb(end)-1);

if projectors
    nk = (k+1)*(k+2)/2;
    PNstarT = zeros(dofPosb(end)-1,nk);
else
    PNstarT = 0;
end

fprintf('Computing local block matrices ...\n')
tic;

rate = zeros(G.cells.num,1);
if ~isempty(src)
    rate(src.cell) = src.rate;
end

for K = 1:nK

    if rem(K,step) == 0
        fprintf('... Calculating local block matrix for cell %d\n', K);
    end 
        
    [AK, bK, dofVec, PNstar] ...
                          = VEM2D_loc(G, K, f, k, alpha, rate(K), mu, rho);

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

[A,b] = VEM2D_bc(G,A,b,bc,k);

end