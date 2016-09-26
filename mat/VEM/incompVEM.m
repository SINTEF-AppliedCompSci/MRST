function [state, varargout] = incompVEM(state, G, S, fluid, varargin)
%   Solves the 2D Poisson equation using a kth order virtual element
%   method.
%
%   SYNOPSIS:
%       [sol, varargout] = VEM2D(G, f, k, bc, varargin)
%
%   DESCRIPTION:
%       Solves the Poisson equation
%
%           -\Delta u = f,
%
%       using the virtual element method of order k. See [1] for details.
%
%   REQUIRED PARAMETERS:
%       G          - 2D MRST grid, with sorted edges, G = sortEdges(G), and
%                    computed VEM geometry, G = computeVEMGeometry(G).
%       f          - Source term. Either a function handle, or a scalar. In
%                    the latter case it is interpreted as a constant
%                    function.
%       k          - Method order. Supported orders are k = 1 and k = 2.
%       bc         - Struct of boundary conditions constructed using
%                    VEM2D_addBC.
%
%   OPTIONAL PARAMETERS:
%       sigma        - G.cells.num x nker matrix of constants for scaling
%                      of the local load terms.
%                      nker = \dim \ker \Pi^\nabla. See [1] for detail.
%       src          - Source term struct constructed using addSource.
%       projectors   - Boolean. If true, matrix representations
%                      of \Pi^\nabla in the monomial basis \mathcal_k(K)
%                      will be added to grid structure G.
%       faceAverages - Boolean. If true, exact face averages of
%                      approximated solution will be calculated
%                      for 1st order VEM. Useful for countour plots.
%       cellAverages - Boolean. If true, exact cell averages of
%                      approximated solution will be calculated
%                      for 1st order VEM. Useful for countour plots.
%
%   RETURNS:
%       sol          - Solution struct. Contans the fileds
%                           * nodeValues, values at the nodes.
%                           * edgeValues, values at the edge
%                             midpoints. Empty for k = 1.
%                           * cellMoments, the first moment (avearge) over
%                             each cell. Empty for k = 1 unless
%                             cellAverages = true.
%
%   OPTIONAL RETURN VALUE:
%       G            - If projectors = true or cellAverages = true, qrid
%                      structure with projectors \Pi^\nabla in the
%                      monomial basis \mathcal_k(K).
%
%   EXAMPLE:
%   
%       G    = cartGrid([10,10]);
%       G    = sortEdges(G)
%       G    = computeVEMGeometry(G);
%       bEdg = find(any(G.faces.neighbors == 0,2));
%       f    = @(X) X(:,1).^2 - X(:,2).^2;
%       bc   = VEM2D_addBC([], boundaryEdges, 'pressure', 0);
%       sol  = VEM2D(G,f,2,bc);
%
%   REFERENCES:
%       [1]     - Ø. S. Klemetsdal: 'The virtual element method as a common
%                 framework for finite element and finite difference
%                 methods - Numerical and theoretical analysis'. MA thesis.
%                 Norwegian University of Science and Technology.
%-----------------------------------------------------------------ØSK-2016-

%{
   Copyright (C) 2016 Øystein Strengehagen Klemetsdal. See COPYRIGHT.txt
   for details.
%}

%%  MERGE INPUT PARAMETRES                                               %%

opt = struct('bc'             , []       , ...
             'src'            , []       , ...
             'srcFunc'        , []       , ...
             'sigma'          , 1        , ...
             'cartGridQ'      , false    , ...
             'faceProjectors' , false    , ...
             'cellProjectors' , false    , ...
             'facePressure'   , false    , ...
             'cellPressure'   , false    , ...
             'linSolve'       , @mldivide, ...
             'matrixOutput'   , false         );
             
opt = merge_options(opt, varargin{:});

%%  CHECK CORRECTNESS OF INPUT                                           %%

assert(G.griddim == 2 || G.griddim == 3, 'Physical dimensin must be 2 or 3.');

k  = S.order;
nN = G.nodes.num;
nF = G.faces.num;
nP = G.cells.num;
% N  = G.nodes.num + G.faces.num*(k-1) + G.cells.num*k*(k-1)/2;

if G.griddim == 2
    nE = 0;
    nk = (k+1)*(k+2)/2;
    NP = diff(G.cells.nodePos) + diff(G.cells.facePos)*(k-1) + k*(k-1)/2;
    N = G.nodes.num + G.faces.num*polyDim(k-2, 1) + G.cells.num*polyDim(k-2, 2);
else
    nE = G.edges.num;
    nk = (k+1)*(k+2)*(k+3)/6;
    NP = diff(G.cells.nodePos) + diff(G.cells.edgePos)*(k-1) ...
       + diff(G.cells.facePos)*k*(k-1)/2 + k*(k^2-1)/6;
    N = G.nodes.num + G.edges.num*polyDim(k-2, 1) + G.faces.num*polyDim(k-2, 2) + G.cells.num*polyDim(k-2, 3);
end
nker = sum(NP - nk);

if isempty(opt.srcFunc)
    opt.srcFunc = 0;
end

if ~isa(opt.srcFunc, 'function_handle')
    assert(numel(opt.srcFunc) == 1 || 0, ...
    'Source function ''srcFunc'' must either be scalar or function handle')
end

assert(any(numel(opt.sigma) == [sum(nker),1]), ...
    'Number of elements in parameter matrix sigma must be 1 or sum(nker)')

if opt.cellPressure
    opt.cellProjectors = true;
end

[mu, rho] = fluid.properties();

%%  ASSEMBLE GLOBAL MATRIX AND COMPUTE RIGHT-HAND SIDE                   %%

[A, rhs] = glob(G, S, opt.src, k, N, mu);

if isempty(opt.bc)
  if k == 1
      n = G.cells.nodes(G.cells.nodePos(1):G.cells.nodePos(2)-1);
      i = n(1);
  else
      if G.griddim == 2  
          i = G.nodes.num + G.faces.num + 1;
      else
          i = G.nodes.num + G.edges.num + G.faces.num + 1;
      end
  end
    A(i,i) = 2*A(i,i); 
else
    [A, rhs] = imposeBC(G, S, opt.bc, k, N, mu, A, rhs);
end

%%  SOLVE LINEAR SYSTEM                                                  %%

% [ii, jj, A] = find(A);
% u = ii > jj; d = ii == jj;
% A = sparse(ii(u), jj(u) , A(u), N, N) + sparse(jj(u), ii(u), A(u), N,N) + sparse(ii(d), jj(d), A(d), N,N); 

p = opt.linSolve(A, rhs);

grav = gravity();
if G.griddim == 3
    if k == 1
        pot = rho*grav(3)*G.nodes.coords(:,3);
    else
        pot = rho*grav(3)*[G.nodes.coords(:,3)   ; G.edges.centroids(:,3); ...
                           G.faces.centroids(:,3); G.cells.centroids(:,3)];
    end
    p = p-pot;
end

%%  UPDATE STATE                                                         %%

state.nodePressure = ...
              full( p( 1:nN)                                             );
state.edgePressure = ...
              full( p((1:nE*(k-1))       + nN)                           );
state.facePressure = ...
              full( p((1:nF*k*(k-1)/2)   + nN + nE*(k-1))                );
state.cellPressure = ...
              full( p((1:nP*k*(k^2-1)/6) + nN + nE*(k-1) + nF*k*(k-1)/2) );

if opt.facePressure && k == 1
    state.facePressure = calculateFacePressure(G, state);
end

if opt.cellPressure && k == 1
    if G.girddim == 2
        state.pressure = calculateCellPressure2D(G,state);
    else
        state.pressure = calculateCellPressure3D(G, state);
    end
end


if opt.matrixOutput
    state.A = A;
    state.rhs = rhs;
end

end

%--------------------------------------------------------------------------

function [A, rhs] = glob(G, S, src, k, N, mu)
    
    P = sparse(1:numel(S.dofVec), S.dofVec, 1, numel(S.dofVec), N);
    A = P'*S.A*P;

    if ~isempty(src)
        if k == 1

            rhs = zeros(G.cells.num,1);
            rhs(src.cell) = src.rate;
            rhs = rldecode(rhs, diff(G.cells.nodePos), 1);
            PiNstar = squeezeBlockDiag(S.PiNstar, diff(G.cells.nodePos), ...
                                       polyDim(1, G.griddim), sum(diff(G.cells.nodePos)))';
            rhs = rhs.*PiNstar(:,1);            
            rhs = mu*P'*rhs;

        else
            rhs = zeros(N,1);
            if G.griddim == 2
                ii = src.cell + G.nodes.num + G.faces.num;
            else
                ii = src.cell + G.nodes.num + G.edges.num + G.faces.num;
            end
            rhs(ii) = mu*src.rate';
        end
        
    else
        rhs = zeros(N,1);
    end
end

%--------------------------------------------------------------------------

function [A, rhs] = imposeBC(G, S, bc, k, N, mu, A, rhs)
    
    nfn = diff(G.faces.nodePos);
    
    if isfield(bc, 'func')
        
        isNeu = find(strcmp(bc.type, 'flux'));
        for i = 1:numel(isNeu)
            f = bc.face{isNeu(i)};
            g = bc.func{isNeu(i)};
            n = G.faces.nodes(mcolon(G.faces.nodePos(f), G.faces.nodePos(f+1)-1));
            if G.griddim == 2
                if k == 1
                    v = bsxfun(@times, bsxfun(@plus, ...
                                reshape(g(G.nodes.coords(n,:)), 2, [])'/6, ...
                                g(G.faces.centroids(f,:))/3), G.faces.areas(f));
                    v = reshape(v', [], 1);
                    rhs = rhs + sparse(n, 1:numel(n), 1, N, numel(n))*v;
                else
                    vn = bsxfun(@times, g(G.nodes.coords(n,:))/6, ...
                                rldecode(G.faces.areas(f), nfn(f), 1));
                    vn = reshape(vn', [], 1);
                    vf = g(G.faces.centroids(f,:))*2/3.*G.faces.areas(f);
                    
                    v = sparse(n, 1:numel(n), 1, N, numel(n))*vn;
                    v(f + G.nodes.num) = vf;
                    rhs = rhs + v;
                end
            else
                
            end

        end
        
        isDir = find(strcmp(bc.type, 'pressure'));
                    
        I = speye(N);
        for i = 1:numel(isDir)
            f = bc.face{isDir(i)};
            g = bc.func{isDir(i)};
            n = G.faces.nodes(mcolon(G.faces.nodePos(f), G.faces.nodePos(f+1)-1));
            n = unique(n);
            if k == 1
                dofVec = n;
                rhs(dofVec) = g(G.nodes.coords(n,:));

            else
                if G.griddim == 2
                    dofVec = [n; f + G.nodes.num];
                    rhs(dofVec) = [g(G.nodes.coords(n,:)); g(G.faces.centroids(f,:))];
                else
                    e = G.faces.edges(mcolon(G.faces.edgePos(f), G.faces.edgePos(f+1)-1));
                    dofVec = [n; e + G.nodes.num; f + G.nodes.num + G.edges.num];
                    rhs(dofVec) = [g(G.nodes.coords(n,:)); ...
                                   g(G.edges.centroids(e,:)); ...
                                   polygonInt3D(G, f, g, 3)./G.faces.areas(f)];
                end
            end
            A(dofVec,:) = I(dofVec,:);
        end
        
    else
       
        NF = diff(G.faces.nodePos);
        
        if G.griddim == 2
            
            isNeu = strcmp(bc.type, 'flux');
            
            if nnz(isNeu)>0
                
                f = bc.face(isNeu);
                v = bc.value(isNeu);
                n = G.faces.nodes(mcolon(G.faces.nodePos(f), G.faces.nodePos(f+1)-1));
                P = sparse(n, 1:numel(n), 1, N, numel(n));
                vn = mu*rldecode(v, NF(f), 1);
                
                if k == 1
                    rhs = rhs + 1/2*P*vn;
                else
                    rhs = rhs + 1/6*P*vn;
                    rhs(f + G.nodes.num) = 2/3*mu*v;
                end
            end
            
            isDir= strcmp(bc.type, 'pressure');
            
            if nnz(isDir)>0
                
                f = bc.face(isDir);
                v = bc.value(isDir);

                n = G.faces.nodes(mcolon(G.faces.nodePos(f), G.faces.nodePos(f+1)-1));
                rhs(n) = rldecode(v, NF(f), 1);

                if k == 1
                    dofVec = n;
                else
                    rhs(f + G.nodes.num) = v;
                    dofVec = [n; f + G.nodes.num];
                end

                I = speye(N);
                A(dofVec,:) = I(dofVec,:);
                
            end

        else
            
        nfe = diff(G.faces.edgePos);
            
        isNeu = strcmp(bc.type, 'flux');
        f = bc.face(isNeu);
        n = G.faces.nodes(mcolon(G.faces.nodePos(f), G.faces.nodePos(f+1)-1));
        v = bc.value(isNeu);
        
        if k == 1
            
            PiNFstar = squeezeBlockDiag(S.PiNFstar, ...
                         NF(G.cells.faces(:,1)), polyDim(k, G.griddim-1), ...
                         sum(NF(G.cells.faces(:,1))));
            pos = [0; cumsum(NF)] + 1;
            v = mu*rldecode(v, NF(f), 1).*PiNFstar(1,mcolon(pos(f), pos(f+1)-1))';
            
            NFf = NF(f);
            ii = rldecode((1:numel(f))', NFf, 1);
            dof = [0; cumsum(NFf(1:end-1))] + 1;
            iiN = mcolon(dof, dof + nfn(f) - 1);
            fDof(iiN) = n;
            v = sum(sparse(fDof, ii, v, N, numel(f)),2);            
            rhs = rhs + v;
            
        else
            rhs(f + G.nodes.num + G.edges.num) = ...
                rhs(f + G.nodes.num + G.edges.num) + mu*v;
        end
        
        

        isDir = strcmp(bc.type, 'pressure');
        f = bc.face(isDir);
        n = G.faces.nodes(mcolon(G.faces.nodePos(f), G.faces.nodePos(f+1)-1));
        v = bc.value(isDir);

        rhs(n) = rldecode(v, nfn(f), 1);

        if k == 1
            dofVec = n;
        else
            e = G.faces.edges(mcolon(G.faces.edgePos(f), G.faces.edgePos(f+1)-1));
            rhs(e + G.nodes.num) = rldecode(v, nfe(f), 1);%./ne;
            rhs(f + G.nodes.num + G.edges.num) = v;
            dofVec = [n; e + G.nodes.num; f + G.nodes.num + G.edges.num];
        end

        I = spdiags(ones(N,1),0,N,N);
        A(dofVec,:) = I(dofVec, :);
        end
            
    end
end
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------

function nk = polyDim(k, dim)

    if k < 0
        nk = 0;
    else
        nk = nchoosek(k+dim,k);
    end
end
