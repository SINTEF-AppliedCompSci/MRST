function [A, b, PNstarT] = VEM2D_glob(G, rock, fluid, k, bc, src, ...
                                     srcFunc, sigma, cartGridQ, projectors)
%   Assmebles the global stiffness matrix and load term for the virtual
%   element method for the 2D Poisson equation.
%
%   SYNOPSIS:
%       [A, b, PNstarT] = VEM2D_glob(G, f, k, bc, sigma, projectors, src)
%
%   DESCRIPTION:
%       Assmebles the global stiffness matrix and load term for the virtual
%       element method of order k on grid G for the 2D Poisson equation 
%
%           -\Delta u = f.
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
%
%   RETURNS:
%       A          - Global stiffness matrix.
%       b          - Global load term.
%       PNstarT    - Projection operators \Pi^\nabla in the monomial basis
%                    \mathcal{M]_k(K). Stored as transpose. See [1] for
%                    details.
%
%   REFERENCES:
%       [1]        - The virtual element method as a common framework for
%                    finite element and finite difference methods -
%                    Numerical and theoretical analysis.
%-----------------------------------------------------------------ØSK-2016-

%{
   Copyright (C) 2016 Øystein Strengehagen Klemetsdal. See Copyright.txt
   for details.
%}

%%  RETRIEVE MONOMIALS, SET FUNCTION SPACE DIMS, CONSTRUCT MAPPINGS      %%

[m, grad_m, int_m] = retrieveMonomials(2,k);

nN = G.nodes.num;
nE = G.faces.num;
nP = G.cells.num;

nk   = (k+1)*(k+2)/2;
NP   = diff(G.cells.nodePos) + diff(G.cells.facePos)*(k-1) + k*(k-1)/2;
nker = NP - nk;
N  = nN + nE*(k-1) + nP*k*(k-1)/2;

dofPosA = [1, cumsum((diff(G.cells.nodePos') + ...
                      diff(G.cells.facePos')*(k-1) + k*(k-1)/2).^2) + 1];
dofPosb = [1, cumsum( diff(G.cells.nodePos') + ...
                      diff(G.cells.facePos')*(k-1) + k*(k-1)/2) + 1];
if numel(sigma) > 1
    sigmaPos = [1, cumsum(nker') + 1];
end

%%  ADJUST INPUT AND OUTPUT PARAMETRES                                   %%

rate = zeros(G.cells.num,1);
if ~isempty(src)
    rate(src.cell) = src.rate;
end

if projectors
    nk = (k+1)*(k+2)/2;
    PNstarT = zeros(dofPosb(end)-1,nk);
else
    PNstarT = 0;
end

K = permTensor(rock,2);
[mu,rho] = fluid.properties();

%%  ASSEMBLE GLOBAL STIFFNESS MATRIX AND LOAD VECTOR                     %%

step = floor(nP/10);

iiA = zeros(1,dofPosA(end)-1);
jjA = zeros(1,dofPosA(end)-1);
AVec = zeros(1,dofPosA(end)-1);

iib = zeros(1,dofPosb(end)-1);
bVec = zeros(1,dofPosb(end)-1);


fprintf('Computing local block matrices ...\n')
tic;

for P = 1:nP

    if rem(P,step) == 0
        fprintf('... Calculating local block matrix for cell %d\n', P);
    end 
    
    if numel(sigma) > 1
        sigmaP = sigma(sigmaPos(P):sigmaPos(P+1)-1);
    else
        sigmaP = sigma;
    end
    
    [AK, bK, dofVec, PNstar] = VEM2D_loc(G, P, K(P,:), mu, rho, k, ...
                    rate(P), srcFunc, sigmaP, cartGridQ, m, grad_m, int_m);

    NP = numel(dofVec);
    
    iiP = repmat(dofVec', NP, 1);
    jjP = repmat(dofVec , NP, 1);
    jjP = jjP(:);
    
    iiA(dofPosA(P):dofPosA(P+1)-1) = iiP;
    jjA(dofPosA(P):dofPosA(P+1)-1) = jjP;
    AVec(dofPosA(P):dofPosA(P+1)-1)= AK(:);

    iib(dofPosb(P):dofPosb(P+1)-1) = dofVec;
    bVec(dofPosb(P):dofPosb(P+1)-1) = bK;
    
    if projectors
        PNstarT(dofPosb(P):dofPosb(P+1)-1,:) = PNstar';
    end
    
end

stop = toc;
fprintf('Done in %f seconds.\n\n', stop);

A = sparse(iiA, jjA, AVec, N, N);
b = sparse(iib, ones(1, numel(iib)), bVec, N, 1);

%%  APPLY BOUNDARY CONDITIONS                                            %%

[A,b] = VEM2D_bc(G,A,b,bc,k);

end