function [A, b, G] = incompVEM3D_glob(G, rock, fluid, k, bc, src, ...
                                srcFunc, sigma, cartGridQ, cellProjectors);
%   Assmebles the global stiffness matrix and load term for the virtual
%   element method for the 2D Poisson equation.
%
%   SYNOPSIS:
%       [A, b, PNstarT] = VEM2D_glob(G, f, k, bc, alpha, projectors, src)
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
%       [1] - Ø. S. Klemetsdal: 'The virtual element method as a common
%             framework for finite element and finite difference methods -
%             Numerical and theoretical analysis'. MA thesis. Norwegian
%             University of Science and Technology.
%-----------------------------------------------------------------ØSK-2016-

%{
   Copyright (C) 2016 Øystein Strengehagen Klemetsdal. See COPYRIGHT.txt
   for details.
%}

%%  RETRIEVE MONOMIALS, SET FUNCTION SPACE DIMS                          %%

[m, ~, ~] = retrieveMonomials(3,k);

nN = G.nodes.num;
nE = G.edges.num;
nF = G.faces.num;
nP = G.cells.num;

nk   = (k+1)*(k+2)*(k+3)/6;
NP   = diff(G.cells.nodePos) + diff(G.cells.edgePos)*(k-1) + ...
       diff(G.cells.facePos)*k*(k-1)/2 + k*(k^2-1)/6;
nker = NP - nk;
N  = nN + nE*(k-1) + nF*k*(k-1)/2 + nP*k*(k^2-1)/6;

%%  COMPUTE FACE PROJECTORS, MONOMIAL DOFS AND SOURCE TERM INTEGRALS     %%


fprintf('Preprocessing ...\n')
tic;

fprintf('... computing face projectors\n');

% G = VEM3D_faceProjectors(G,k);

if k == 2    
    fprintf('... computing source term integrals\n');
    
    if isa(f, 'function_handle')
        IFf = polygonInt3D(G,1:G.faces.num,f,k+1);
        ICf = polyhedronInt(G,1:G.cells.num,f,k+1);
    else
        IFf = f*G.faces.areas;
        ICf = f*G.cells.volumes;
    end
    G.faces.('fFaceIntegrals') = IFf;
    G.cells.('fCellIntegrals') = ICf;
end

stop = toc;
fprintf('Done in %f seconds.\n\n', stop);


%%  CREATE MAPPINGS                                                      %%

dofPosA = [1, cumsum((NP.^2)') + 1];
dofPosb = [1, cumsum( NP'    ) + 1];

%%  ADJUST INPUT AND OUTPUT PARAMETRES                                   %%

if numel(sigma) > 1
    sigmaPos = [1, cumsum(nker') + 1];
end

rate = zeros(G.cells.num,1);
if ~isempty(src)
    rate(src.cell) = src.rate;
end

if cellProjectors
    PNstarT = zeros(dofPosb(end)-1,nk);
else
    PNstarT = 0;
end

K = permTensor(rock,2);
[mu, rho] = fluid.properties();

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
    
    [AP, bP, dofVec, PNstar] = incompVEM3D_loc(G, P, K(P,:), mu, rho, k, rate(P), ...
                                            srcFunc, sigmaP, cartGridQ, m);

    NP = numel(dofVec);
    
    iiP = repmat(dofVec', NP, 1);
    jjP = repmat(dofVec , NP, 1);
    jjP = jjP(:);
    
    iiA(dofPosA(P):dofPosA(P+1)-1) = iiP;
    jjA(dofPosA(P):dofPosA(P+1)-1) = jjP;
    AVec(dofPosA(P):dofPosA(P+1)-1)= AP(:);

    iib(dofPosb(P):dofPosb(P+1)-1) = dofVec;
    bVec(dofPosb(P):dofPosb(P+1)-1) = bP;
    
    if cellProjectors
        PNstarT(dofPosb(P):dofPosb(P+1)-1,:) = PNstar';
    end
    
end

stop = toc;
fprintf('Done in %f seconds.\n\n', stop);

A = sparse(iiA, jjA, AVec, N, N);
b = full(sparse(iib, ones(1, numel(iib)), bVec, N, 1));

if cellProjectors
    G.cells.('PNstarT') = PNstarT;
    PNstarPos = [1, cumsum(diff(G.cells.nodePos')            + ...
                           diff(G.cells.edgePos')*(k-1)      + ...
                           diff(G.cells.facePos')*k*(k-1)/2  + ...
                           k*(k^2-1)/6) + 1                  ];
    G.cells.('PNstarPos') = PNstarPos;
end

%%  APPLY BOUNDARY CONDITIONS                                            %%

[A,b] = VEM3D_bc(G,A,b,bc,k);

end