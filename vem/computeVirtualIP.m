function S = computeVirtualIP(G, rock, k, varargin)
%Compute local inner products for the frist- or second order virtual
%element method (VEM).
%
% SYNOPSIS:
%   S = computeVirtualIP(G, rock, k)
%   S = computeVirtualIP(G, rock, k, pn1, 'vn1', ...)
%
% REQUIRED PARAMETERS:
%   G    - Grid structure as described by grid_structure.
%
%   rock - Rock data structure with valid field 'perm'.  The
%          permeability is assumed to be in measured in units of
%          metres squared (m^2).  Use function 'darcy' to convert from
%          (milli)darcies to m^2, e.g.,
%
%          perm = convertFrom(perm, milli*darcy)
%
%          if the permeability is provided in units of millidarcies.
%
%          The field rock.perm may have ONE column for a scalar
%          permeability in each cell, TWO/THREE columns for a diagonal
%          permeability in each cell (in 2/3 D) and THREE/SIX columns
%          for a symmetric full tensor permeability.  In the latter
%          case, each cell gets the permeability tensor
%
%          K_i = [ k1  k2 ]      in two space dimensions
%                [ k2  k3 ]
%
%          K_i = [ k1  k2  k3 ]  in three space dimensions
%                [ k2  k4  k5 ]
%                [ k3  k5  k6 ]
%
%   k    - Method order. A k-th order method vil recover k-th order
%          pressure fields exactly. Supported values are 1 and 2.
%
% OPTIONAL PARAMETERS:
%
%   innerProduct - The choice of stability term in the inner
%                  product. String. Default value = 'ip_simple'.
%                  Supported values are:
%                    - 'ip_simple'  : 'Standard' VEM stability term equal
%                                      to trace(K)h^(dim-2) I.
%
%                    - 'ip_qfamily' : Parametric family of inner products.
% 
%                    - 'ip_fem'     : Inner product resembling the finite
%                                     element mehtod on regular Cartesian
%                                     grids.
%                    - 'ip_fd'      : Inner product resembling a
%                                     combination of two finite difference
%                                     stencils, one using the regular
%                                     Cartesian coordinate axes (Fc), and
%                                     one using the cell diagonals (Fd) as
%                                     axes on regular Cartesian grids.
%
%  sigma          - Extra parameters to inner product
%                   ip_qfamily. Must be either a single scalar value, or
%                   nker values per cell, where nker is the dim of the
%                   nullspace of the projeciton operator \Pi^\nabla.
%
%  w              - Extra parameter to the inner product ip_fd,
%                   corresponding to the weighting of the two FD stencils,
%                   wFc + (1-w)Fd. Positive scalar vale.
%
%  invertBlocks   - Method by which to invert a sequence of
%                   small matrices that arise in the
%                   discretisation. String.
%                   Supported values are:
%                     - MATLAB : Use the MATLAB function mldivide
%                                (backslash) (the default).
%
%                     - MEX    : Use two C-accelerated MEX functions to
%                                extract and invert, respectively, the
%                                blocks along the diagonal of a sparse
%                                matrix.  This method is often faster by a
%                                significant margin, but relies on being
%                                able to build the required MEX functions.
%
%  trans          - Fluxes can alternatively be reconstructed from the
%                   VEM pressure field using the TPFA or MPFA scheme.
%                   String. Supported values are 'mpfa' and 'tpfa'.
%
% RETURNS:
%  S - Pressure linear system structure having the following fields:
%        - A       : Block diagonal matrix with the local inner products
%                    on the diagonal.
%        - order   : Order k of the VEM method.
%        - ip      : Inner product name.
%        - PiNstar : Block diagonal matrix with the Projection operators
%                    \Pi^\nabla in the monomial basis for each cell.
%        - PiNFstar: Block diagonal matrix with the Projection operators
%                    \Pi^\nabla in the monomial basis for each face of
%                    each cell. Empty if G.griddim = 2.
%        - faceCoords: Local 2D coordinate systems for each face if
%                    G.griddim = 3.
%        - dofVec  : Map from local to global degrees of freedom.
%        - T, transType: If MPFA or TPFA scheme is to be used for flux
%                    reconstructiom, T and transType are the corresponding
%                    transmissibilites and type.
%
% SEE ALSO:
%   incompVEM, darcy, permTensor.

%{
Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}

%   Written by Øystein Strengehagen Klemetsdal, SINTEF/NTNU, 2016.

%%  MERGE INPUT PARAMETRES                                               %%

opt = struct('innerProduct', 'ip_simple', ...
             'sigma'       , []         , ...
             'w'           , []         , ...
             'invertBlocks', 'MATLAB'   , ...
             'trans'       , []          );
opt = merge_options(opt, varargin{:});

%%  CFUNCTION SPACE DIMENSIONS

%   Polynomial space dimension
nk = polyDim(k, G.griddim);

ncn = diff(G.cells.nodePos);
ncf = diff(G.cells.facePos);
if G.griddim == 3; nce = diff(G.cells.edgePos);
else nce = zeros(G.cells.num,1); G.edges.num = 0; end

%   Number of dofs for each face
if G.griddim == 3
    nfn = diff(G.faces.nodePos);
    nfe = diff(G.faces.edgePos);
    NF  = nfn ...
        + nfe*polyDim(k-2, G.griddim-2) ...
        +     polyDim(k-2, G.griddim-1);
end

%   Number of dofs for each cell
NP = ncn                           ...
   + nce*polyDim(k-2, G.griddim-2) ...
   + ncf*polyDim(k-2, G.griddim-1) ...
   +     polyDim(k-2, G.griddim);

%   Dimension of \ker \Pi^\nabla
nker = NP - nk;

%   Total number of dofs
N = G.nodes.num                           ...
  + G.edges.num*polyDim(k-2, G.griddim-2) ...
  + G.faces.num*polyDim(k-2, G.griddim-1) ...
  + G.cells.num*polyDim(k-2, G.griddim  );

%%  CHECK INPUT

if any(strcmp(opt.innerProduct, {'ip_fem', 'ip_fd'})) && ...
   ((~any(strcmp(G.type, 'cart')) && k > 1) || G.griddim > 2)
    warning(['innerProduct ''ip_fem'' and ''ip_fd'' only valid for 2D 1st order and', ...
             'Cartesian grid. Using ip_simple instead.']);
    opt.innerProduct = 'ip_simple';
end

ipNames = {'ip_simple' , 'ip_qfamily', 'ip_fem', 'ip_fd'};

assert(any(strcmp(opt.innerProduct, ipNames)), ...
                        'Unknown inner product ''%s''.', opt.innerProduct);

if ~isempty(opt.sigma)
    assert((strcmp(opt.innerProduct, 'ip_qfamily') ...
                                   && numel(opt.sigma) == sum(nker)) || ...
           (any(strcmp(opt.innerProduct, {'ip_fem', 'ip_fd'})) ...
                                   && numel(opt.sigma) == G.cells.num));
end

% if G.griddim == 3 && strcmp(opt.invertBlocks, 'MATLAB')
%     warning(['Computing VEM inner products in 3D without ', ...
%              '''invertBlocks'' set to ''MEX'' can be very slow. ', ...
%              'Press ENTER to continue, or ctrl-C to abort']);
%     pause;
% end

%%  EXTRACT PREMABILITY, FACES, EDGES AND NODES

K = permTensor(rock, G.griddim);

%   Faces for each cell.
f    = G.cells.faces(:,1);
fn   = G.faces.normals(f,:);
fSgn = (-ones(numel(f),1)).^(G.faces.neighbors(f,1) ...
       ~= rldecode((1:G.cells.num)', diff(G.cells.facePos), 1)); 
fn   = bsxfun(@times, fn, fSgn);
if size(f,1) == 1; f = f'; end

if G.griddim == 2

    %   Nodes for each face of each cell.
    n   = G.faces.nodes(mcolon(G.faces.nodePos(f), G.faces.nodePos(f+1)-1));
    if size(n,1) == 1; n = n'; end
    n   = reshape(n,2,[])';
    n(fSgn == -1,:) = n(fSgn == -1,2:-1:1);
    n   = n(:,1);
    
else
    
    %   Edges for each face of each cell
    eNum = mcolon(G.faces.edgePos(f), G.faces.edgePos(f+1)-1);
    e    = G.faces.edges(eNum);
    if size(e,1) == 1; e = e'; end
    en = G.faces.edgeNormals(eNum,:);
    
    %   Nodes for each edge of each face of each cell
    n = G.edges.nodes(mcolon(G.edges.nodePos(e), G.edges.nodePos(e+1)-1));
    if size(n,1) == 1; n = n'; end
    n = reshape(n,2,[])';
    n(G.faces.edgeSign(eNum) == -1,:) = n(G.faces.edgeSign(eNum) == -1, 2:-1:1);
    n = n(:,1);
end

%%  CALCULATE MULTIPOINT TRANSMISSIBILITY FOR FLUX CALCULATIONS          %%

if strcmp(opt.trans, 'mpfa')
    T = computeMultiPointTrans(G, rock);
elseif strcmp(opt.trans, 'tpfa')
    T = computeTrans(G, rock);
else
    T = [];
end

if G.griddim == 2

    %%  CALCULATE LOCAL STIFFNESS MATRICES FOR EACH CELL
    
    %   Coordinates for degrees of freedom.
    if k == 1
        x = G.nodes.coords(n,:);
    else
        x = [G.nodes.coords(n,:); G.faces.centroids(f,:)];
    end
    
    Kmat = reshape(K', 2, [])';
    
    %   Calculate B and D matrices.
    [B, D] = computeBD2D(G.cells.centroids, G.cells.diameters, ...
                         G.cells.volumes, ncn, ncf, fn,        ...
                         G.faces.areas(f), G.cells.facePos, x, ...
                         numel(n), G.cells.nodePos, Kmat, NP, k);
    
    %   Calculate projection operators in monomial (star) and VEM bases.
    M = B*D;
    [ii, jj] = blockDiagIndex(repmat(nk, [G.cells.num ,1]));
    kk = sub2ind(size(M), ii, jj);
    
    if strcmp(opt.invertBlocks, 'MEX')
        PiNstar = sparse(ii, jj, invv(full(M(kk)), repmat(nk, [G.cells.num, 1])))*B;
    else
        PiNstar = M\B;
    end
    PiN = D*PiNstar;

    clear B D;

    SS = stabilityTerm(G, K, PiN, NP, nker, opt);
    
    M(1:nk:end,:) = 0;
    I = speye(size(PiN,1));
    A = PiNstar'*M*PiNstar + (I-PiN)'*SS*(I-PiN);

    %   Make solution struct.
    S = makeSolutionStruct(G, NP, k, A, T, PiNstar, [], [], [], opt);
    
else
    
    %%  CALCULATE PROJECTION OPERATORS FOR EACH FACE
    
    nkf = polyDim(k, G.griddim-1);
    
    %   Compute local coordinates for each face.
    [v1, v2, xf, ecf, enf] = faceCoorSys(G);
     
    %   Project permeability tensors onto each face.
    Kmat = rldecode(K, diff(G.cells.facePos), 1); Kmat = Kmat';
    [ii,jj] = blockDiagIndex(repmat(3,[size(Kmat,2), 1]));
    Kmat = sparse(ii, jj, Kmat(:));
    F = [v1; v2]; F = F(:,f);
    [ii, jj] = blockDiagIndex(3*ones(numel(f),1), 2*ones(numel(f),1));
    F = sparse(ii, jj, F(:));
    Kmat = squeezeBlockDiag(F'*Kmat*F, repmat(2, [numel(f), 1]), 2*numel(f), 2);
      
    %   Coordinates for degrees of freedom. 
    iin = mcolon(G.faces.nodePos(f), G.faces.nodePos(f+1)-1);
    iie = mcolon(G.faces.edgePos(f), G.faces.edgePos(f+1)-1);
    if k == 1; xx = xf(iin,:); else xx = [xf(iin,:); ecf(iie,:)]; end
    
    ePos = diff(G.faces.edgePos); ePos = ePos(f); ePos = [1;cumsum(ePos)+1];
    nPos = diff(G.faces.nodePos); nPos = nPos(f); nPos = [1;cumsum(nPos)+1];
    
    %   Calculate B and D matrices.
    [BF, DF] = computeBD2D(zeros(numel(f),2), G.faces.diameters(f), ...
                           G.faces.areas(f), nfn(f), nfe(f), ...
                           enf(iie,:), G.edges.lengths(e), ePos, xx, ...
                           numel(n), nPos, Kmat, NF(f), k);

    %   Calculate projection operators in monomial (star) and VEM bases.
    MF = BF*DF;
    [ii, jj] = blockDiagIndex(repmat(nkf, [numel(f) ,1]));
    kk = sub2ind(size(MF), ii, jj);
    PiNFstar = sparse(ii, jj, invv(full(MF(kk)), ...
                      repmat(nkf, [numel(f), 1]))   )*BF;
%     PiNFstar = MF\BF;
    
    clear BF DF;
    
    %%  CALCULATE PROJECTION OPERATORS AND LOCAL STIFFNESS MATRICES

    [ii, jj] = blockDiagIndex(ones(numel(f),1), nfe(f));
    Ie = sparse(ii, jj, 1);
    
    [ii, jj] = blockDiagIndex(ones(G.cells.num,1), ncf);
    If = sparse(ii, jj, 1);
    
    %%  CALCULATE D MATRICES

    %   Fix face normal signs and normalize.
    fSgn = (-ones(numel(f),1)).^(G.faces.neighbors(f,1) ...
              ~= rldecode((1:G.cells.num)', diff(G.cells.facePos), 1)); 
    fn = bsxfun(@times, G.faces.normals(f,:), fSgn./G.faces.areas(f));
    
    %   Express monomial coordinates in face coordinates.
    fc  = rldecode(G.faces.centroids(f,:), nfn(f), 1);
 
    [jj, ii] = blockDiagIndex(3*ones(numel(f),1), nfn(f));
    x = (G.nodes.coords(n,:)-fc)'     ; x = sparse(ii, jj, x(:));
    x   = squeezeBlockDiag(x*F, nfn(f), sum(nfn(f)), 2);
    ec  = (G.edges.centroids(e,:)-fc)'; ec  = sparse(ii, jj, ec(:));
    ec  = squeezeBlockDiag(ec*F, nfn(f), sum(nfn(f)), 2);
    en  = en'                         ; en  = sparse(ii, jj, en(:));
    en  = squeezeBlockDiag(en*F, nfn(f), sum(nfn(f)), 2);
    enx = en(:,1).*G.edges.lengths(e);

    if k == 1
        
        %   For k = 1, the entries of D are the values of the monomials at
        %   the veritces of each cell.
        
        [alpha, beta, gamma] = retrieveMonomials(k, G.griddim);
        xMon = bsxfun(@rdivide, G.nodes.coords(G.cells.nodes,:) ...
                               - rldecode(G.cells.centroids, ncn,1), ...
                                 rldecode(G.cells.diameters, ncn, 1));
        mVals = (bsxfun(@power, xMon(:,1), alpha) ...
              .*bsxfun(@power, xMon(:,2), beta ) ...
              .*bsxfun(@power, xMon(:,3), gamma))';
          
        [jj, ii] = blockDiagIndex(nk*ones(G.cells.num,1), NP);
        D = sparse(ii, jj, mVals(:));
                
    else
        
        %   For k = 2 the entries of D are the values of the monomials at
        %   the vertices and edge centroids of each cells, the first
        %   moments of each monomial over each face of each cell, and the
        %   first moments of each monomial over each cell.
        
        [alpha, beta, gamma] = retrieveMonomials(k, G.griddim);
        
        %   Evaluate monomials at cell vetrices.
        xMon = bsxfun(@rdivide, G.nodes.coords(G.cells.nodes,:) ...
                              - rldecode(G.cells.centroids, ncn, 1), ...
                                rldecode(G.cells.diameters, ncn, 1));
                            
        mn =   bsxfun(@power, repmat(xMon(:,1),1, nk), alpha) ...
             .*bsxfun(@power, repmat(xMon(:,2),1, nk), beta ) ...
             .*bsxfun(@power, repmat(xMon(:,3),1, nk), gamma );
         
        %   Evaluate monomials at cell edge centroids.
        ecMon = bsxfun(@rdivide, G.edges.centroids(G.cells.edges,:) ...
                              - rldecode(G.cells.centroids, nce, 1), ...
                                rldecode(G.cells.diameters, nce, 1));
        
        me =   bsxfun(@power, repmat(ecMon(:,1),1, nk), alpha) ...
             .*bsxfun(@power, repmat(ecMon(:,2),1, nk), beta ) ...
             .*bsxfun(@power, repmat(ecMon(:,3),1, nk), gamma );
        
        %   Evaluate first moment of each linear monomial over each cell
        %   face using centroid rule. (First moment of each linear monomial
        %   over each cell is zero).
        
        %   Use centroid rule for the linear monomials.
        fcMon = bsxfun(@rdivide, G.faces.centroids(f,:) ...
                              - rldecode(G.cells.centroids, ncf, 1), ...
                                rldecode(G.cells.diameters, ncf, 1));
         
        linfInt = bsxfun(@power, repmat(fcMon(:,1),1,3), alpha(2:4) ) ...
                 .*bsxfun(@power, repmat(fcMon(:,2),1, 3), beta (2:4) ) ...
                 .*bsxfun(@power, repmat(fcMon(:,3),1, 3), gamma(2:4) );
            
        %   Use divergence theorem twice for the bilinear and quadratic
        %   monomials.
             
        %   Repeat cell centroids in number of faces per cell.
        ccf = rldecode(G.cells.centroids, ncf, 1);

        %   Express x-x_P, y-y_P, and z-z_P in face coordinates.
        cx = trinomialExpansion(v1(1,f)',v2(1,f)', ...
                                     G.faces.centroids(f,1) - ccf(:,1), 1);
        cy = trinomialExpansion(v1(2,f)',v2(2,f)', ...
                                     G.faces.centroids(f,2) - ccf(:,2), 1);
        cz = trinomialExpansion(v1(3,f)',v2(3,f)', ...
                                     G.faces.centroids(f,3) - ccf(:,3), 1);
        
        %   Express bilinear monomials in local face coordinates.
        alpha = [1 0 0]; beta  = [0 1 0];
        [alphaBi, betaBi, c6] = polyProducts(cx, cy, alpha, beta);
        [~      , ~     , c7] = polyProducts(cx, cz, alpha, beta);
        [~      , ~     , c9] = polyProducts(cy, cz, alpha, beta);
        cBi = [c6', c7', c9'];
        
        %   Divide by x-exponent.
        alphaBi = alphaBi+1;
        cBi = bsxfun(@rdivide, cBi, alphaBi');
        
        %   Express quadratic monimials in local face coordinates.
        c5  = trinomialExpansion(v1(1,f)', v2(1,f)', ...
                                     G.faces.centroids(f,1) - ccf(:,1), 2);
        c8  = trinomialExpansion(v1(2,f)', v2(2,f)', ...
                                     G.faces.centroids(f,2) - ccf(:,2), 2);
        c10 = trinomialExpansion(v1(3,f)', v2(3,f)', ...
                                     G.faces.centroids(f,3) - ccf(:,3), 2);
        cQuad = [c5', c8', c10'];
        
        %   Divide by x-exponent.
        alphaQuad = [2 1 1 0 0 0];
        betaQuad  = [0 1 0 2 1 0];
        alphaQuad = alphaQuad + 1;
        cQuad = bsxfun(@rdivide, cQuad, alphaQuad');

        %   Integrate each monimial over each face of each cell using
        %   three-point Gauss-Lobatto on each edge.
        
        pos = [1;cumsum(nfn(f))+1];
        iim  = 1:size(x,1); jjm = iim;
        jjm(1:end-1) = jjm(2:end);
        jjm(cumsum(pos(2:end)-pos(1:end-1))) ...
                           = iim([1;cumsum(pos(2:end-1) - pos(1:end-2))+1]);
        
        mVals = bsxfun(@power, repmat([x(:,1); ec(:,1)],1,6), alphaQuad)...
              .*bsxfun(@power, repmat([x(:,2); ec(:,2)],1,6), betaQuad);
        mVals = bsxfun(@times, (mVals(iim,:) + mVals(jjm,:))/6 ...
                              + mVals(size(x,1)+1:end,:)*2/3, enx)';
        [jj, ii] = blockDiagIndex(size(alphaQuad,2)*ones(numel(f),1), nfn(f));
        mVals = sparse(ii, jj, mVals(:));
        
        quadfInt = Ie*mVals*reshape(cQuad, size(mVals,2), 3);
        
        mVals = bsxfun(@power, repmat([x(:,1); ec(:,1)],1,3*3), alphaBi)...
              .*bsxfun(@power, repmat([x(:,2); ec(:,2)],1,3*3), betaBi);
        mVals = bsxfun(@times, (mVals(iim,:) + mVals(jjm,:))/6 ...
                              + mVals(size(x,1)+1:end,:)*2/3, enx)';
        [jj, ii] = blockDiagIndex(size(alphaBi,2)*ones(numel(f),1), nfn(f));
        mVals = sparse(ii, jj, mVals(:));
        
        bifInt = Ie*mVals*reshape(cBi, size(mVals,2), 3);
        
        cd = rldecode(G.cells.diameters, diff(G.cells.facePos), 1);
        fInt = bsxfun(@rdivide, ...
                      [quadfInt(:,1), bifInt(:,1:2), ...
                       quadfInt(:,2), bifInt(:,3) , quadfInt(:,3)],...
                       G.faces.areas(f).*cd.^2);
        
        %   Calculate first moments of each monomial over each cell.
        
        %   Express antiderivatives of bilinear monomials in face
        %   coordinates.
        cy = [cy, zeros(numel(f), polyDim(k, G.griddim) - polyDim(k-1,G.griddim))];
        cz = [cz, zeros(numel(f), polyDim(k, G.griddim) - polyDim(k-1,G.griddim))];        
        c5 = [zeros(numel(f), polyDim(k-1, G.griddim)-1), c5];
        c8 = [zeros(numel(f), polyDim(k-1, G.griddim)-1), c8];
        
        alpha = [1 0 0 2 1 1 0 0 0];
        beta  = [0 1 0 0 1 0 2 1 0];
        [alphaBi, betaBi, c6] = polyProducts(c5, cy, alpha, beta);
        [~      , ~     , c7] = polyProducts(c5, cz, alpha, beta);
        [~      , ~     , c9] = polyProducts(c8, cz, alpha, beta);
        
        alpha = [0 1 0 2 1 0 3 2 1 0];
        beta  = [0 0 1 0 1 2 0 1 2 3];
        I = bsxfun(@eq, alphaBi', repmat(alpha, 9*9, 1)) & bsxfun(@eq, betaBi', repmat(beta, 9*9, 1));
        c6 = c6*I; c7 = c7*I; c9 = c9*I;
        alphaBi  = alpha; betaBi = beta;
        
        c6 = bsxfun(@times, c6, fn(:,1)/2);
        c7 = bsxfun(@times, c7, fn(:,1)/2);
        c9 = bsxfun(@times, c9, fn(:,2)/2);
         
        alphaBi = alphaBi +1;
        cBi = bsxfun(@rdivide, [c6', c7', c9'], alphaBi');
        
        %   Express antiderivatives of quadratic monimials in face
        %   coordinates.
        c5  = trinomialExpansion(v1(1,f)', v2(1,f)', ...
                                     G.faces.centroids(f,1) - ccf(:,1), 3);
        c8  = trinomialExpansion(v1(2,f)', v2(2,f)', ...
                                     G.faces.centroids(f,2) - ccf(:,2), 3);
        c10 = trinomialExpansion(v1(3,f)', v2(3,f)', ...
                                     G.faces.centroids(f,3) - ccf(:,3), 3);
        
        c5  = bsxfun(@times, c5 , fn(:,1)/3);
        c8  = bsxfun(@times, c8 , fn(:,2)/3);
        c10 = bsxfun(@times, c10, fn(:,3)/3);
        
        alphaQuad = [3 2 2 1 1 1 0 0 0 0];
        betaQuad  = [0 1 0 2 0 1 3 2 1 0];
        alphaQuad = alphaQuad + 1;
        cQuad = bsxfun(@rdivide, [c5', c8', c10'], alphaQuad');
        
        %   Integrate each monomial over each cell using four-point
        %   Gauss-Lobatto on each edge.
        
        %   Compute Gauss-lobatto quadrature points and express in local
        %   coordinates.
        eVec = G.nodes.coords(G.edges.nodes(2:2:end),:)...
             - G.nodes.coords(G.edges.nodes(1:2:end),:);
        xq1 = G.edges.centroids - .5*sqrt(1/5)*eVec;
        xq2 = G.edges.centroids + .5*sqrt(1/5)*eVec;
        
        [jj, ii] = blockDiagIndex(3*ones(numel(f),1), nfn(f));
        
        xq1 = (xq1(e,:)-fc)'     ; xq1 = sparse(ii, jj, xq1(:));
        xq1 = squeezeBlockDiag(xq1*F, nfn(f), sum(nfn(f)), 2);
        
        xq2 = (xq2(e,:)-fc)'     ; xq2 = sparse(ii, jj, xq2(:));
        xq2 = squeezeBlockDiag(xq2*F, nfn(f), sum(nfn(f)), 2);
        
        %   Integrate each monomial over each edge of each face of each
        %   cell.
        
        nn = size(x,1);
        mVals = bsxfun(@power, ...
                   repmat([x(:,1); xq1(:,1); xq2(:,1)],1,nk), alphaQuad)...
              .*bsxfun(@power, ...
                   repmat([x(:,2); xq1(:,2); xq2(:,2)],1,nk), betaQuad );
        mVals = bsxfun(@times, (mVals(iim,:) + mVals(jjm,:))/12 ...
               + (mVals(nn + 1:2*nn,:) + mVals(2*nn+1:3*nn,:))*5/12, enx)';
        [jj, ii] = blockDiagIndex(size(alphaQuad,2)*ones(numel(f),1), nfn(f));
        mVals = sparse(ii, jj, mVals(:));

        quadcInt = If*Ie*mVals*reshape(cQuad, size(mVals,2), 3);
        
        mVals = bsxfun(@power, ...
                    repmat([x(:,1); xq1(:,1); xq2(:,1)],1,numel(alphaBi)), alphaBi)...
              .*bsxfun(@power, ...
                    repmat([x(:,2); xq1(:,2); xq2(:,2)],1,numel(betaBi)), betaBi );

        mVals = bsxfun(@times, (mVals(iim,:) + mVals(jjm,:))/12 ...
               + (mVals(nn + 1:2*nn,:) + mVals(2*nn+1:3*nn,:))*5/12, enx)';
        [jj, ii] = blockDiagIndex(size(alphaBi,2)*ones(numel(f),1), nfn(f));
        mVals = sparse(ii, jj, mVals(:));
        
        bicInt = If*Ie*mVals*reshape(cBi, size(mVals,2), 3);
        
        cInt = bsxfun(@rdivide, ...
            [quadcInt(:,1), bicInt(:,1:2), quadcInt(:,2), bicInt(:,3), quadcInt(:,3)], ...
            G.cells.volumes.*G.cells.diameters.^2);
    
        %   Assmeble D matrices.
        dof = [0; cumsum(NP(1:end-1))] + 1;
        iiN = mcolon(dof, dof + ncn - 1);
        iiE = mcolon(dof + ncn, dof + ncn + nce - 1);
        iiF = mcolon(dof + ncn + nce, dof + ncn + nce + ncf - 1);
        iiP = mcolon(dof + ncn + nce + ncf, dof + ncn + nce + ncf);

        D([iiN, iiE, iiF, iiP], :) ...
            = [mn; me; ones(numel(f), 1), linfInt, fInt; ...
               ones(G.cells.num,1), zeros(G.cells.num, 3), cInt];
        D = D';
        [jj, ii] = blockDiagIndex(nk*ones(G.cells.num,1), NP);
        D = sparse(ii, jj, D(:));

    end
    
    %% CALCULATE B MATRICES
    
    %   Calculate integrals of \nabla \phi^i K \nabla m^\alpha over each
    %   cell using
    %   \int_P \nabla \phi^i K \nabla m^\alpha dx = 
    %         \int_\partial \Pi^\nabla \phi^i K \nabla m^\alpha \cdot n ds.
    
    %   Calculate K n_f
    Kmat  = reshape(K', 3, [])';
    [jj, ii] = blockDiagIndex(3*ones(G.cells.num,1), ncf);
    fn = fn'; c = sparse(ii,jj,fn(:))*Kmat;
    
    if k == 1
        
        %   All integrands are linear, and we use centroid rule to evaluate
        %   integrals.
        
        c = rldecode(c, NF(f), 1);
        PiNFs = squeezeBlockDiag(PiNFstar, NF(f), polyDim(k, 2), sum(NF(f)));
        c = bsxfun(@times, c, PiNFs(1,:)'.*rldecode(G.faces.areas(f), NF(f), 1));
        
        %   Map to global coordinates
        ii = rldecode((1:numel(f))', NF(f), 1);
        int2 = sparse(ii, n, c(:,1));
        int3 = sparse(ii, n, c(:,2));
        int4 = sparse(ii, n, c(:,3));
        
        %   For each cell, sum all face integrals
        int2 = (If*int2)'; int2 = int2(:);
        int3 = (If*int3)'; int3 = int3(:);
        int4 = (If*int4)'; int4 = int4(:);
        int  = [int2(:), int3(:), int4(:)];

        vec  = repmat(G.nodes.num,G.cells.num,1);
        vec  = [0; cumsum(vec(1:end-1))];
        cDof = G.cells.nodes + rldecode(vec, NP,1);
        int  = int(cDof,:);
        
        %   Assemble B matrices.
        BT = zeros(sum(NP), nk);
        cd = rldecode(G.cells.diameters, NP, 1);
        BT(:,2:end) = bsxfun(@rdivide, int, cd);
        BT(:,1)     = rldecode(1./NP, NP, 1);
        B = BT';
        [ii, jj] = blockDiagIndex(nk*ones(G.cells.num,1), NP);
        B = sparse(ii, jj, B(:));

    else

        c2 = c(:,1); c3 = c(:,2); c4 = c(:,3);
        c5 = c(:,1)*2; c8 = c(:,2)*2; c10 = c(:,3)*2;

        %   Express x-x_P, y-y_P and z-z_P in face coordinates.
        cx = trinomialExpansion(v1(1,f)', v2(1,f)', ...
                                       G.faces.centroids(f,1)-ccf(:,1), 1);
        cy = trinomialExpansion(v1(2,f)', v2(2,f)', ...
                                       G.faces.centroids(f,2)-ccf(:,2), 1);
        cz = trinomialExpansion(v1(3,f)', v2(3,f)', ...
                                       G.faces.centroids(f,3)-ccf(:,3), 1);

        zer = zeros(numel(f),3);
        c5  = bsxfun(@times, cx(:,[3,1,2]), c5 ); c5  = [c5 , zer];
        c6  = bsxfun(@times, cy(:,[3,1,2]), c2 ) ...
            + bsxfun(@times, cx(:,[3,1,2]), c3 ); c6  = [c6 , zer];
        c7  = bsxfun(@times, cz(:,[3,1,2]), c2 ) ...
            + bsxfun(@times, cx(:,[3,1,2]), c4 ); c7  = [c7 , zer];
        c8  = bsxfun(@times, cy(:,[3,1,2]), c8 ); c8  = [c8 , zer];
        c9  = bsxfun(@times, cz(:,[3,1,2]), c3 ) ...
            + bsxfun(@times, cy(:,[3,1,2]), c4 ); c9  = [c9 , zer];
        c10 = bsxfun(@times, cz(:,[3,1,2]), c10); c10 = [c10, zer];

        %   Put coefficients in a suitable format.
        PiNFs = squeezeBlockDiag(PiNFstar', NF(f), sum(NF(f)), polyDim(k,2));
        c5  = rldecode(c5 , NF(f), 1);    
        c6  = rldecode(c6 , NF(f), 1);
        c7  = rldecode(c7 , NF(f), 1);
        c8  = rldecode(c8 , NF(f), 1);
        c9  = rldecode(c9 , NF(f), 1);
        c10 = rldecode(c10, NF(f), 1);

        %   Multiply coeffiecients of the monomials
        a = [0 1 0 2 1 0];
        b = [0 0 1 0 1 2];
        [alpha510, beta510, c5 ] = polyProducts(c5 , PiNFs, a, b);
        [~       , ~      , c6 ] = polyProducts(c6 , PiNFs, a, b);
        [~       , ~      , c7 ] = polyProducts(c7 , PiNFs, a, b);
        [~       , ~      , c8 ] = polyProducts(c8 , PiNFs, a, b);
        [~       , ~      , c9 ] = polyProducts(c9 , PiNFs, a, b);
        [~       , ~      , c10] = polyProducts(c10, PiNFs, a, b);
        
        alpha = [0 1 0 2 1 0]; beta = [0 0 1 0 1 2];
        alpha = repmat(alpha, 1,6); beta = repmat(beta,1 ,6);
        fd = bsxfun(@power, repmat(rldecode(G.faces.diameters(f), ...
               NF(f), 1), 1, polyDim(k, G.griddim-1)^2), alpha + beta);
        
        c5 = c5./fd; c6 = c6./fd; c7 = c7./fd; c8 = c8./fd; c9 = c9./fd; c10 = c10./fd;
        
        alpha = [0 1 0 2 1 0 3 2 1 0 4 3 2 1 0];
        beta  = [0 0 1 0 1 2 0 1 2 3 0 1 2 3 4];
        I = bsxfun(@eq, alpha510', repmat(alpha, 6*6, 1)) & bsxfun(@eq, beta510', repmat(beta, 6*6, 1));
        c5 = c5*I; c6 = c6*I; c7 = c7*I; c8 = c8*I; c9 = c9*I; c10 = c10*I;
        alpha510  = alpha; beta510 = beta;
       
        %   Add 1 to all x-coordinate exponents, and divide by result.
        alpha510 = alpha510 + 1;

        c5  = bsxfun(@rdivide, c5 , alpha510)';
        c6  = bsxfun(@rdivide, c6 , alpha510)';
        c7  = bsxfun(@rdivide, c7 , alpha510)';
        c8  = bsxfun(@rdivide, c8 , alpha510)';
        c9  = bsxfun(@rdivide, c9 , alpha510)';
        c10 = bsxfun(@rdivide, c10, alpha510)';
        
        %   Put in block diagonal format.
        
        [ii, jj] = blockDiagIndex(numel(alpha510)*ones(numel(f),1),  NF(f));
        c5  = sparse(ii, jj, c5 (:));
        c6  = sparse(ii, jj, c6 (:));
        c7  = sparse(ii, jj, c7 (:));
        c8  = sparse(ii, jj, c8 (:));
        c9  = sparse(ii, jj, c9 (:));
        c10 = sparse(ii, jj, c10(:));

        %   Integrate monomials over each edge of each face using four-point
        %   (m^5-m^10) Gauss-Lobatto quadratures.

        m5m10 = bsxfun(@power, repmat([x(:,1); xq1(:,1); xq2(:,1)],1,15), alpha510)...
               .*bsxfun(@power, repmat([x(:,2); xq1(:,2); xq2(:,2)],1,15), beta510);

        pos = [1;cumsum(nfn(f))+1];
        ii = 1:size(x,1); jj = ii;
        jj(1:end-1) = jj(2:end);
        jj(cumsum(pos(2:end)-pos(1:end-1))) = ii([1;cumsum(pos(2:end-1) - pos(1:end-2))+1]);

        m5m10 = (Ie*bsxfun(@times, (m5m10(ii,:) + m5m10(jj,:))/12 ...
                    + (m5m10(nn + 1:2*nn,:) + m5m10(2*nn+1:3*nn,:))*5/12, enx))';
        [jj, ii] = blockDiagIndex(15*ones(numel(f),1), ones(numel(f),1));
        m5m10 = sparse(ii, jj, m5m10(:));

        %   Map to global dofs and sum for each cell
        int5  = squeezeBlockDiag(m5m10*c5 , NF(f), 1, sum(NF(f)));   
        int6  = squeezeBlockDiag(m5m10*c6 , NF(f), 1, sum(NF(f)));
        int7  = squeezeBlockDiag(m5m10*c7 , NF(f), 1, sum(NF(f)));
        int8  = squeezeBlockDiag(m5m10*c8 , NF(f), 1, sum(NF(f)));
        int9  = squeezeBlockDiag(m5m10*c9 , NF(f), 1, sum(NF(f)));
        int10 = squeezeBlockDiag(m5m10*c10, NF(f), 1, sum(NF(f)));

        NFf = NF(f);
        ii = rldecode((1:numel(f))', NFf, 1);
        dof = [0; cumsum(NFf(1:end-1))] + 1;
        iiN = mcolon(dof, dof + nfn(f) - 1);
        iiE = mcolon(dof + nfn(f), dof + nfn(f) + nfe(f) - 1);
        iiF = mcolon(dof + nfn(f) + nfe(f), dof + nfn(f) + nfe(f));

        int24 = zeros(sum(NF(f)),3);
        int24(cumsum(NF(f)),:) = bsxfun(@times, c, G.faces.areas(f));

        int2 = int24(:,1); int3 = int24(:,2); int4 = int24(:,3);

        fDof([iiN, iiE, iiF]) = [n; ...
                                 e + G.nodes.num; ...
                                 f + G.nodes.num + G.edges.num];

        int2  = sparse(ii, fDof, int2 , numel(f), N);
        int3  = sparse(ii, fDof, int3 , numel(f), N);
        int4  = sparse(ii, fDof, int4 , numel(f), N);
        int5  = sparse(ii, fDof, int5 , numel(f), N);
        int6  = sparse(ii, fDof, int6 , numel(f), N);
        int7  = sparse(ii, fDof, int7 , numel(f), N);
        int8  = sparse(ii, fDof, int8 , numel(f), N);
        int9  = sparse(ii, fDof, int9 , numel(f), N);
        int10 = sparse(ii, fDof, int10, numel(f), N);
        
        int2  = (If*int2)' ; int2  = int2 (:);
        int3  = (If*int3)' ; int3  = int3 (:);
        int4  = (If*int4)' ; int4  = int4 (:);
        int5  = (If*int5)' ; int5  = int5 (:);
        int6  = (If*int6)' ; int6  = int6 (:);
        int7  = (If*int7)' ; int7  = int7 (:);
        int8  = (If*int8)' ; int8  = int8 (:);
        int9  = (If*int9)' ; int9  = int9 (:);
        int10 = (If*int10)'; int10 = int10(:);

        %   Pick out the dofs for each cell
        dof = [0; cumsum(NP(1:end-1))]+1;
        iiN = mcolon(dof, dof + ncn -1)';
        iiE = mcolon(dof + ncn, dof + ncn + nce -1)';
        iiF = mcolon(dof + ncn + nce, dof + ncn + nce + ncf - 1)';
        cDof = zeros(sum(NP), 1);
        cDof([iiN; iiE; iiF]) = [G.cells.nodes; ...
                                 G.cells.edges + G.nodes.num; ...
                                 f             + G.nodes.num + G.edges.num];
        cDof = cDof(cDof~= 0,:);
        vec = repmat(N,G.cells.num,1);
        vec = [0; cumsum(vec(1:end-1))];
        cDof = cDof + rldecode(vec, NP-1,1);

        int = [int2, int3, int4, int5, int6, int7, int8, int9, int10];
        int = int(cDof,:);

        %   Build B matrices
        BT = zeros(sum(NP), nk);
        vec = [0; cumsum(NP)] + 1;
        ii = mcolon(vec(1:end-1), vec(2:end)-2);
        cdi = rldecode(G.cells.diameters, NP-1, 1);
        BT(ii,2:end) = [bsxfun(@rdivide, int(:, 1:3), cdi), ...  
                        bsxfun(@rdivide, int(:, 4:9), cdi.^2)];



        vec = zeros(G.cells.num,nk);
        vec(:, [1,5:nk]) = [ones(G.cells.num,1), ...
                           bsxfun(@times, -2*[K(:,1:3), K(:,5:6), K(:,9)], ...
                           G.cells.volumes./G.cells.diameters.^2)];
        BT(cumsum(NP),:) = BT(cumsum(NP),:) + vec;
        B = BT';
        
        [ii, jj] = blockDiagIndex(nk*ones(G.cells.num,1), NP);
        B = sparse(ii, jj, B(:));

    end
    
    M = B*D;
    
    [ii, jj] = blockDiagIndex(repmat(nk, [G.cells.num ,1]));
    kk = sub2ind(size(M), ii, jj);
    
    if strcmp(opt.invertBlocks, 'MEX')
        PiNstar = sparse(ii, jj, invv(full(M(kk)), repmat(nk, [G.cells.num, 1])))*B;
    else
        PiNstar = M\B;
    end
    PiN = D*PiNstar; 
    
    SS = stabilityTerm(G, K, PiN, NP, nker, opt);
    M(1:nk:end,:) = 0;
    I = speye(size(PiN,1));
    A = PiNstar'*M*PiNstar + (I-PiN)'*SS*(I-PiN);

    %   Make solution struct.
    S = makeSolutionStruct(G, NP, k, A, T, PiNstar, PiNFstar, v1, v2, opt);
    
end

end
    
%--------------------------------------------------------------------------

function [B, D] = computeBD2D(cc, cd, cv, ncn, ncf, fn, fa, fPos, x, nn, nPos, K, NP, k)
%   Calculate B and D matrices in 2D.

nk = polyDim(k, 2);
nc = numel(cv);
nf = numel(fa);

%%  CALCULATE D MATRICES

%   D_{(i, \alpha)} is the ith degree of freedom of monomial m^\alpha.

if k == 1
    
    %   For k = 1, D is the monomial values at the cell nodes.
    
    alpha = [0 1 0]; beta = [0 0 1];
    xMon = bsxfun(@rdivide, x - rldecode(cc, ncn, 1), rldecode(cd, ncn,1));
    D = (bsxfun(@power, xMon(:,1), alpha).*bsxfun(@power, xMon(:,2), beta))';
    [jj, ii] = blockDiagIndex(nk*ones(nc,1), NP);
    D = sparse(ii, jj, D(:));
    
else
    
    %   For k = 2, D is the monomial values at the cell nodes and face
    %   centroids, and the first moment over the cell.
    
    %   Evaluate anti-derivative of monimals.
    alpha = [0 1 0 2 1 0]; beta = [0 0 1 0 1 2]; alpha = alpha + 1;
    xMon = bsxfun(@rdivide, x - repmat(rldecode(cc, ncn, 1),2,1), ...
                                repmat(rldecode(cd, ncn,1),2,1));
    mVals = bsxfun(@rdivide, bsxfun(@power, xMon(:,1), alpha)...
                           .*bsxfun(@power, xMon(:,2), beta ), alpha);
    
    %   Integrate monomial over each face using three-point Gauss-Lobatto.
    ii = 1:nn; jj = ii;
    jj(1:end-1) = jj(2:end);
    jj(cumsum(ncn)) = ii([1;cumsum(ncn(1:end-1))+1]);
    intD = bsxfun(@times, ( mVals(ii,:) + mVals(jj,:) )/6 ...
                          + mVals(nn+1:end,:)*2/3, ...
                            fn(:,1)                   );

    %   Sum for each face of each cell.
    [ii, jj] = blockDiagIndex(ones(nc,1), ncf);
    If = sparse(ii, jj, 1);
    intD = bsxfun(@times, If*intD, cd./cv);
    
    %   Assmeble matrices.
    nodeDof = mcolon([1;cumsum(NP(1:end-1))+1],[1;cumsum(NP(1:end-1))+1]+ncn-1);
    edgeDof = nodeDof + rldecode(ncn, ncn, 1)';
    D = zeros(sum(NP), nk);
    
    alpha = [0 1 0 2 1 0]; beta = [0 0 1 0 1 2];
    
    D([nodeDof, edgeDof],:) ...
       = bsxfun(@power, xMon(:,1), alpha).*bsxfun(@power, xMon(:,2), beta);
    D(cumsum(NP),:) = intD;
    D = D';
    [jj, ii] = blockDiagIndex(nk*ones(nc,1), NP);
    D = sparse(ii, jj, D(:));
    
end

%%  CALCULATE B MATRICES

%   B_{(\alpha, i)} is the integral
%       a^P(\phi^i, m^\alpha)
%           =  \int_\partial \phi^i K \nabla m^\alpha \cdot n ds
%             -\int_P \phi^i \nabla K \nabla m^\alpha  dx.

if k == 1
    
    %   Only linear monomials, use centroide rule.
    
    %   Calculate K n_f
    [ii,jj] = blockDiagIndex(2*ones(size(cc,1),1));
    K = K'; K = sparse(ii,jj,K(:));
    
    fn = bsxfun(@rdivide, fn, rldecode(cd, ncf, 1));
    [jj, ii] = blockDiagIndex(2*ones(size(cc,1),1), ncf);
    fn = fn'; fn = sparse(ii,jj,fn(:));

    %   Use centroid rule to evaluate integrals.
    BT = .5*fn*K;
    ii = 1:nf; jj = ii;
    jj(2:end) = jj(1:end-1);
    jj([1;cumsum(fPos(2:end-1)    ...
                -fPos(1:end-2))+1]) ...
     = ii(cumsum(fPos(2:end)-fPos(1:end-1)));
    BT = BT(ii,:) + BT(jj,:);
 
    %   Assemble matrices (first row is the integral of each basis function
    %   over the boundary.
    BT = [.5*(fa(ii) + fa(jj)), squeezeBlockDiag(BT,ncf, nf, nk-1)];
    B = BT';
    [ii, jj] = blockDiagIndex(nk*ones(nc,1), NP);
    B = sparse(ii, jj, B(:));
    
else
    
    %   Evaluate monomial gradients.
    alphax = [0 0 1 0 0]; beta  = [0 1 0 1 2]; cx = [1 0 2 1 0];
    alpha  = [1 0 2 1 0]; betay = [0 0 0 0 1]; cy = [0 1 0 1 2];
    gm = bsxfun(@times, ...
         [bsxfun(@power, xMon(:,1), alphax).*bsxfun(@power, xMon(:,2), beta), ...
          bsxfun(@power, xMon(:,1), alpha) .*bsxfun(@power, xMon(:,2), betay)], [cx, cy]);
    
    %   Put in sparse block format.
    ii = repmat((1:5*sum(ncn + ncf))',2,1);
    jj = repmat(1:2,5*sum(ncn + ncf),1);
    add = repmat([rldecode((0:2:2*(nc-1))', ncn,1); ...
                  rldecode((0:2:2*(nc-1))', ncf,1)], 5,1);
    jj = bsxfun(@plus,jj,add); jj = jj(:);
    intB = sparse(ii, jj, gm(:), 5*sum(ncn+ ncf), 2*nc)*K;
    
    %   Dot product by length-weighted face normals.
    ii = 1:nn; jj = ii;
    jj(2:end) = jj(1:end-1);
    jj(nPos(1:end-1)) = ii(nPos(2:end)-1);
    
    iin = repmat(1:nn, 1, 5) + rldecode(0:(nn+nf):4*(nn+nf), nn*ones(1,5), 2);
    iif =  iin + nn;
    
    intB = sum(intB([iin, iin, iif],:).*...
               [repmat(fn,5,1); repmat(fn(jj,:),5,1); repmat(fn,5,1)], 2);
    
    %   Evaluate line integrals using three-point Gauss-Lobatto.    
    intB = [reshape((intB(1:numel(iin)) + intB(numel(iin)+1:2*numel(iin)))/6, nn, 5);
            reshape(intB(2*numel(iin)+1:end)*2/3, nn, 5)];
    intB = bsxfun(@rdivide, intB, repmat(rldecode(cd,diff(nPos),1),2,1));

    %   Assmble matrices.
    BT = zeros(sum(NP), nk);
    BT([nodeDof, edgeDof],2:nk) = intB;
    K = reshape(K', 4, [])';
    vec = zeros(nc,6);
    vec(:, [1,4:6]) = [ones(nc,1), ...
                       bsxfun(@times, -2*[K(:,1),K(:,2), K(:,4)], ...
                       cv./cd.^2)];
    BT(cumsum(NP),:) = BT(cumsum(NP), :) + vec;
    
    B = BT';
    [ii, jj] = blockDiagIndex(nk*ones(nc,1), NP);
    B = sparse(ii, jj, B(:));

end

end

%--------------------------------------------------------------------------

function [v1, v2, x, ec, en] = faceCoorSys(G)
%   Construct local coordinate system for each face. For face f, the
%   mapping from global to face coordinates is
%   (x-G.faces.centroids(f,:))*[v1(:,f), v2(:,f)].
%   x, ec and en are node coordinates, edge centroids and edge normals in
%   local face coordinates
    
    %   Face edges.
    e  = G.faces.edges;   
    en = G.faces.edgeNormals;
    en = bsxfun(@times, en, G.edges.lengths(e));

    %   Face nodes.
    n   = G.edges.nodes(mcolon(G.edges.nodePos(e),G.edges.nodePos(e+1)-1));
    n   = reshape(n,2,[])';
    n(G.faces.edgeSign == -1,:) = n(G.faces.edgeSign == -1,2:-1:1);
    n   = n(:,1);
    
    %   Node coordinates.
    x = G.nodes.coords(n,:);
    
    %   Create mapping.
    v1 = (x(G.faces.nodePos(1:end-1)+1,:) - x(G.faces.nodePos(1:end-1),:));
    v1 = bsxfun(@rdivide, v1, sqrt(sum(v1.^2,2)));
    v2 = cross(G.faces.normals,v1,2);
    v2 = bsxfun(@rdivide, v2, sqrt(sum(v2.^2,2)));
    v1 = v1'; v2 = v2';
    T(:,[1:2:2*G.faces.num, 2:2:2*G.faces.num]) = [v1, v2];
    [ii, jj] = blockDiagIndex(3*ones(G.faces.num,1), 2*ones(G.faces.num,1));
    T = sparse(ii, jj, T(:));
    
    nfn = diff(G.faces.nodePos);
    nfe = diff(G.faces.edgePos);
    
    fc = rldecode(G.faces.centroids, nfn, 1);
    
    %   Map from global to local cooridnates.
    [jj, ii] = blockDiagIndex(3*ones(G.faces.num,1), nfn);
    x = (x-fc)'                       ; x = sparse(ii, jj, x(:));
    x   = squeezeBlockDiag(x*T, nfn, sum(nfn), 2);
    ec  = (G.edges.centroids(e,:)-fc)'; ec  = sparse(ii, jj, ec(:));
    ec  = squeezeBlockDiag(ec*T, nfn, sum(nfn), 2);
    en  = en'                         ; en  = sparse(ii, jj, en(:));
    en  = squeezeBlockDiag(en*T, nfn, sum(nfn), 2);


end

%--------------------------------------------------------------------------

function SS = stabilityTerm(G, K, PiN, NP, nker, opt)
%   Construct stability term. 
%
%   ip_simple : Identity matrix multiplied by trace of permeability tensor
%               and cell diameter to the power of dim-2.
%
%   ip_qfamily: Stability term Q \Sigma Q^T, where Q is an orthogonal basis
%               for \ker \Pi^\nabla, and \Sigma is an n_\ker x n_\ker
%               diagonal matrix, with diagonal entries specified by input
%               sigma.
%   ip_fem    : FEM for Cartesian grids. Only supported in 2D 1st order.
%   ip_fd     : FD for Cartesian grids. Only supported in 2D 1st order.

    if G.griddim == 2; iiK = [1 4]; else iiK = [1 5 9]; end
    
    innerProduct = opt.innerProduct;
    
    switch innerProduct

        case 'ip_simple'
            
           
            SS = spdiags(rldecode(G.cells.diameters.^(G.griddim-2)  ...
                                  .*sum(K(:,iiK),2)/numel(iiK),NP,1), ...
                                  0, sum(NP), sum(NP)              );
        
        case 'ip_qfamily'
            
            Q = zeros(sum(nker.*NP),1);
            PiNPos = [1; cumsum(NP.^2) + 1];
            QPos   = [1; cumsum(NP.*nker)+1];
            ii = blockDiagIndex(NP);
            PiNvec = full(PiN(ii));

            for P = 1:G.cells.num 
                QP = null(reshape(PiNvec(PiNPos(P):PiNPos(P+1)-1), ...
                          NP(P), NP(P))                               );
                Q(QPos(P):QPos(P+1)-1) = QP(:);
            end

            [ii,jj] = blockDiagIndex(NP, nker);
            Q = sparse(ii, jj, Q, sum(NP), sum(nker));
            
            sigma = opt.sigma;
            if isempty(sigma)
                sigma = rldecode(G.cells.diameters.^(G.griddim-2) ...
                                 .*sum(K(:,iiK),2),nker,1);
            end

            SS = Q*spdiags(sigma, 0, sum(nker), sum(nker))*Q';
            
        case {'ip_fem', 'ip_fd'}
            
            xx = G.nodes.coords(G.cells.nodes,:);
            [ii, jj] = blockDiagIndex(diff(G.cells.nodePos), ones(G.cells.num,1));
            x = reshape(xx(:,1), 4, []);
            y = reshape(xx(:,2), 4, []);
            
            hx = abs(max(x,[], 1)- min(x,[],1))'/2;
            hy = abs(max(y,[], 1)- min(y,[],1))'/2;

            Q = sqrt(9./(4*rldecode(hx.*hy, 4*ones(G.cells.num,1),1)))...
                                     .*repmat([-1,1,-1,1]', G.cells.num,1);
            Q = sparse(ii,jj,Q);
            
%             if isempty(opt.sigma)
%                 sigma = ones(G.cells.num,1);
%             end
            
            factor = 1; w = opt.w;
            if strcmp(innerProduct, 'ip_fd')
                if isempty(w)
                    w = 1;
                end
                factor = 3*w;
            end
            Lambda = spdiags(factor*3*(K(:,1)./hx.^2 + K(:,4)./hy.^2), ...
                                              0, G.cells.num, G.cells.num);
            SS = Q*((Q'*Q)\Lambda/(Q'*Q))*Q';
    end
    
end

%--------------------------------------------------------------------------

function S = makeSolutionStruct(G, NP, k, A, T, PiNstar, PiNFstar, v1, v2, opt)
%   Make solution struct.

    ncn = diff(G.cells.nodePos);
    ncf = diff(G.cells.facePos);
    vec = [1; cumsum(NP(1:end-1)) + 1];
    iiN = mcolon(vec, vec + ncn-1);
        
    if G.griddim == 2
        
            %   Faces for each cell.
        f    = G.cells.faces(:,1);
        fn   = G.faces.normals(f,:);
        fSgn = (-ones(numel(f),1)).^(G.faces.neighbors(f,1) ...
               ~= rldecode((1:G.cells.num)', diff(G.cells.facePos), 1)); 
        if size(f,1) == 1; f = f'; end


        %   Nodes for each face of each cell.
        n   = G.faces.nodes(mcolon(G.faces.nodePos(f), G.faces.nodePos(f+1)-1));
        if size(n,1) == 1; n = n'; end
        n   = reshape(n,2,[])';
        n(fSgn == -1,:) = n(fSgn == -1,2:-1:1);
        n   = n(:,1)';

                
        iiF = mcolon(vec + ncn, vec + ncn + ncf*polyDim(k-2, 1) -1);
        iiP = mcolon(vec + ncn + ncf*polyDim(k-2, 1), ...
                     vec + ncn + ncf*polyDim(k-2, 1) + polyDim(k-2, 2) -1);
        if k == 1
            dofVec([iiN, iiF, iiP]) = n;
        else
            dofVec([iiN, iiF, iiP]) ...
                = [n,...
                   G.cells.faces(:,1)' + G.nodes.num, ...
                   (1:G.cells.num) + G.nodes.num + G.faces.num*polyDim(k-2, 2)];
        end
        
    else
        
        nce = diff(G.cells.edgePos);
        iiE = mcolon(vec + ncn, vec + ncn + nce*polyDim(k-2, 1) - 1);
        iiF = mcolon(vec + ncn + nce*polyDim(k-2, 1), ...
                     vec + ncn + nce*polyDim(k-2, 1) + ncf*polyDim(k-2, 2) - 1);
        iiP = mcolon(vec + ncn + nce*polyDim(k-2, 1) + ncf*polyDim(k-2, 2), ...
                     vec + ncn + nce*polyDim(k-2, 1) + ncf*polyDim(k-2, 2) + polyDim(k-2,3)-1);

        if k == 1
            dofVec([iiN, iiE, iiF, iiP]) = G.cells.nodes';
        else
            dofVec([iiN, iiE, iiF, iiP]) ...
                = [G.cells.nodes', ...
                   G.cells.edges' + G.nodes.num, ...
                   G.cells.faces(:,1)' + G.nodes.num ...
                   + G.edges.num*polyDim(k-2, 1), ...
                   (1:G.cells.num) + G.nodes.num ...
                   + G.edges.num*polyDim(k-2, 1) + G.faces.num*polyDim(k-2, 2)];
        end

    end
    
    S.A          = A;
    S.ip         = opt.innerProduct;
    S.dofVec     = dofVec;
    S.PiNstar    = PiNstar;
    S.PiNFstar   = PiNFstar;
    S.faceCoords = [v1(:), v2(:)];
    S.order      = k;
    S.T          = T;
    S.transType  = opt.trans;
    
end
    
%--------------------------------------------------------------------------

function coeff = trinomialExpansion(a, b, c, n)
%   Expands polynomial on the form (ax + by + c)^n to the form
%   \sum_i coeff(i)*x^alpha(i)*y^beta(i)

    if n == 0
        alpha = 0; beta = 0; gamma = 0;
    elseif n == 1
        alpha = [1,0,0]; beta = [0,1,0]; gamma = [0,0,1];
    elseif n == 2
        alpha = [2,1,1,0,0,0]; beta = [0,1,0,2,1,0]; gamma = [0,0,1,0,1,2];
    else
        alpha = [3 2 2 1 1 1 0 0 0 0];
        beta  = [0 1 0 2 0 1 3 2 1 0];
        gamma = [0 0 1 0 2 1 0 1 2 3];
    end
    
    r = size(a,1);     
    coeff = repmat(factorial(n)...
            ./(factorial(alpha).*factorial(beta).*factorial(gamma)), r, 1);
    coeff = coeff.*bsxfun(@power, a,repmat(alpha,r,1))...
         .*bsxfun(@power, b,repmat(beta,r,1))...
         .*bsxfun(@power, c,repmat(gamma,r,1));
    
end

%--------------------------------------------------------------------------

function [alpha, beta, coeff] = polyProducts(coeff1,coeff2,alph, bet)
%   Calculates the product of polynomials on the form (\sum_i
%   coeff1(i)*x^alph(i)*y^bet(i))*(sum_i coeff2(i)*x^alph(i)*y^bet(i)) into
%   the form \sum_i coeff(i)*x^alpha(i)*y^beta(i).

    [r,c] = size(coeff1);
    cPos  = 1:c:c*c+1;
    coeff = zeros(r, cPos(end)-1);
    alpha = zeros(1, cPos(end)-1);
    beta  = zeros(1, cPos(end)-1);
    for i = 1:c
        coeff(:, cPos(i):cPos(i+1)-1) = coeff1(:, [i:end, 1:i-1]).*coeff2;
        alpha(cPos(i):cPos(i+1)-1) = alph + alph([i:end, 1:i-1]);
        beta(cPos(i):cPos(i+1)-1) = bet + bet([i:end, 1:i-1]);
    end
end

function [alpha, beta, gamma] = retrieveMonomials(k, dim)
    if dim == 2
        gamma = [];
        if k == 1
            alpha = [0 1 0];
            beta  = [0 0 1];
        else
            alpha = [0 1 0 2 1 0];
            beta  = [0 0 1 0 1 2];
        end
    else
        if k == 1
            alpha = [0 1 0 0];
            beta  = [0 0 1 0];
            gamma = [0 0 0 1];
        else
            alpha = [0 1 0 0 2 1 1 0 0 0];
            beta  = [0 0 1 0 0 1 0 2 1 0];
            gamma = [0 0 0 1 0 0 1 0 1 2];
        end
    end
end
           
%--------------------------------------------------------------------------