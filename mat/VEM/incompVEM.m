function state = incompVEM(state, G, S, fluid, varargin)
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
             'facePressure'   , false    , ...
             'cellPressure'   , false    , ...
             'conservativeFlux', false, ...
             'linSolve'       , @mldivide, ...
             'matrixOutput'   , false         );
             
opt = merge_options(opt, varargin{:});

%%  CHECK CORRECTNESS OF INPUT                                           %%

assert(G.griddim == 2 || G.griddim == 3, 'Physical dimension must be 2 or 3.');

if G.griddim == 2
    G.edges.num = 0;
end

N = G.nodes.num + G.edges.num*polyDim(S.order - 2, G.griddim - 2) ...
                + G.faces.num*polyDim(S.order - 2, G.griddim - 1) ...
                + G.cells.num*polyDim(S.order - 2, G.griddim    );

[mu, rho] = fluid.properties();

%%  ASSEMBLE GLOBAL MATRIX AND COMPUTE RIGHT-HAND SIDE                   %%

[A, rhs] = assembleSystem(state, G, S, fluid, opt);

%%  SOLVE LINEAR SYSTEM                                                  %%

p = opt.linSolve(A, rhs);

gvec = gravity();
if G.griddim == 3
    if S.order == 1
        pot = rho*gvec(3)*G.nodes.coords(:,3);
    else
        pot = rho*gvec(3)*[G.nodes.coords(:,3)   ; G.edges.centroids(:,3); ...
                           G.faces.centroids(:,3); G.cells.centroids(:,3)];
    end
    p = p-pot;
end

%%  UPDATE PRESSURE STATE                                                %%

state.nodePressure ...
               = full(p(1:G.nodes.num));
state.edgePressure ...
               = full(p((1:G.edges.num*polyDim(S.order-2, G.griddim-2)) ...
                           + G.nodes.num));
state.facePressure ...
               = full(p((1:G.faces.num*polyDim(S.order-2, G.griddim-1)) ...
                           + G.nodes.num ...
                           + G.edges.num*polyDim(S.order-2, G.griddim-2)));
state.pressure     ...
               = full(p((1:G.cells.num*polyDim(S.order-2, G.griddim  )) ...
                           + G.nodes.num ...
                           + G.edges.num*polyDim(S.order-2, G.griddim-2)...
                           + G.faces.num*polyDim(S.order-2, G.griddim-1)));
                       
if opt.facePressure
    state.facePressure = facePressure(state, G, S, 1:G.faces.num);
end

if (opt.cellPressure || ~isempty(S.T)) && isempty(state.pressure)
    state.pressure = cellPressure(state, G, S);
    if isempty(state.facePressure)
        state.facePressure(boundaryFaces(G)) = facePressure(state, G, S, boundaryFaces(G));
    end
end
    
if ~isempty(S.T)
    state.flux = computeFlux(state, G, S, opt.bc, mu)/mu;
    if opt.conservativeFlux
        state.flux = conserveFlux(state, G, opt.bc);
    end
end

if opt.matrixOutput
    state.A   = A;
    state.rhs = rhs;
end

end

function [A, rhs] = assembleSystem(state, G, S, fluid, opt)
    
if G.griddim == 2
    G.edges.num = 0;
end

[mu, rho] = fluid.properties();

N = G.nodes.num + G.edges.num*polyDim(S.order - 2, G.griddim - 2) ...
                + G.faces.num*polyDim(S.order - 2, G.griddim - 1) ...
                + G.cells.num*polyDim(S.order - 2, G.griddim    );

    [A, rhs] = glob(G, S, opt.src, N, mu);

%     totmob = dynamic_quantities(state, fluid);
%     
%     tm = zeros(N,1);
%     
%     ii = rldecode((1:G.nodes.num)', diff(G.nodes.cellPos),1);
%     jj = G.nodes.cells;
%     P = sparse(ii, jj,1);
%     tm(1:G.nodes.num) = P*totmob./sum(P,2);
%     
%     tm((1:G.edges.num*polyDim(k-2,G.griddim-2)) + G.nodes.num) ...
%                   = sparse(G.cells.edges, 1:numel(G.cells.edges), 1)*totmob;
%               
%     tm((1:G.faces.num*polyDim(k-2, G.griddim-1)) + G.nodes.num + G.edges.num ) ...
%                  = sparse(G.cells.faces, 1:numel(G.cells.faces), 1)*totmob;
%              
%     tm((1:G.cells.num*polyDim(k-2, G.griddim-1)) + G.nodes.num + G.edges.num +G.faces.num) ...
%                  = sparse(G.cells.faces, 1:numel(G.cells.faces), 1)*totmob;
%    
%     A = spdiags(tm, 0, N, N)\A;
    
    pressure_bc = ~isempty(opt.bc) && any(strcmpi('pressure', opt.bc.type));

    if ~isempty(opt.bc)
        [A, rhs] = imposeBC(A, rhs, G, S, opt.bc, N, mu);
    end

    if ~pressure_bc
      if S.order == 1
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
    end
    


end

%--------------------------------------------------------------------------

function [A, rhs] = glob(G, S, src, N, mu)
    
    P = sparse(1:numel(S.dofVec), S.dofVec, 1, numel(S.dofVec), N);
    A = P'*S.A*P;

    if ~isempty(src)
        if S.order == 1

            rhs = zeros(G.cells.num,1);
            rhs(src.cell) = src.rate;
            rhs = rldecode(rhs, diff(G.cells.nodePos), 1);
            PiNstar = squeezeBlockDiag(S.PiNstar, diff(G.cells.nodePos), ...
                 polyDim(S.order, G.griddim), sum(diff(G.cells.nodePos)))';
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

function [A, rhs] = imposeBC(A, rhs, G, S, bc, N, mu)
    
    nfn = diff(G.faces.nodePos);
    
    if isfield(bc, 'func')
        
        isNeu = find(strcmp(bc.type, 'flux'));
        for i = 1:numel(isNeu)
            f = bc.face{isNeu(i)};
            g = bc.func{isNeu(i)};
            n = G.faces.nodes(mcolon(G.faces.nodePos(f), G.faces.nodePos(f+1)-1));
            if G.griddim == 2
                if S.order == 1
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
            if S.order == 1
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
                
                if S.order == 1
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

                if S.order == 1
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
        
        if S.order == 1
            
            PiNFstar = squeezeBlockDiag(S.PiNFstar, ...
                         NF(G.cells.faces(:,1)), polyDim(S.order, G.griddim-1), ...
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

        if S.order == 1
            dofVec = n;
        else
            e = G.faces.edges(mcolon(G.faces.edgePos(f), G.faces.edgePos(f+1)-1));
            rhs(e + G.nodes.num) = rldecode(v, nfe(f), 1);%./ne;
            rhs(f + G.nodes.num + G.edges.num) = v;
            dofVec = [n; e + G.nodes.num; f + G.nodes.num + G.edges.num];
        end

        I = speye(N);
        A(dofVec,:) = I(dofVec, :);
        end
            
    end
end

function p = facePressure(state, G, S, f)
     
    if G.griddim == 2
        
        n = G.faces.nodes(mcolon(G.faces.nodePos(f),G.faces.nodePos(f+1)-1));
        [ii, jj] = blockDiagIndex(ones(numel(f),1), 2*ones(numel(f),1));
        p = 1/2*sparse(ii, jj, 1)*state.nodePressure(n);

    else
        
        cf = G.cells.faces(:,1);
        nfn = diff(G.faces.nodePos);
        c = squeezeBlockDiag(S.PiNFstar, nfn(cf),...
              polyDim(S.order, G.griddim-1), sum(cf));
        c = c(1,:)';

        ii = rldecode((1:numel(cf))', nfn(cf), 1);
        jj = 1:sum(nfn(cf));
        I = sparse(ii,jj,1);
                
        e = G.faces.edges(mcolon(G.faces.edgePos(cf), G.faces.edgePos(cf+1)-1));
        n = G.edges.nodes(mcolon(G.edges.nodePos(e),G.edges.nodePos(e+1)-1));
        n = reshape(n, 2, [])';
        n(G.faces.edgeSign(f) == -1,:) = n(G.faces.edgeSign(f) == -1, 2:-1:1);
        n = n(:,1);
        
        p = full(I*(c.*state.nodePressure(n)));
        
        ii = cf;
        jj = 1:numel(cf);
        I = sparse(ii, jj, 1);
        p = I*p./sum(I,2);        
        p = p(f);
        
    end
end

function p = cellPressure(state, G, S)
    
    ncn = diff(G.cells.nodePos);
    c = squeezeBlockDiag(S.PiNstar, ncn, polyDim(S.order, G.griddim), sum(ncn));
    c = c(1,:)';
    ii = rldecode((1:G.cells.num)', ncn, 1);
    jj = 1:numel(G.cells.nodes);
    I = sparse(ii,jj,1);
    p = full(I*(c.*state.nodePressure(G.cells.nodes)));
    
end

%--------------------------------------------------------------------------

function flux = computeFlux(state, G, S, bc, mu)
    
    bf = boundaryFaces(G);
    
    ii = [(1:G.cells.num)'; G.cells.num + bf];
    p = [state.pressure; state.facePressure(bf)];
    flux = cellFlux2faceFlux(G, S.T.T(:, ii) * p);
    
    if ~isfield(bc, 'func')
    
        isNeu = strcmp(bc.type, 'flux');
        f = bc.face(isNeu);
        flux(f) = bc.value(isNeu)*mu;

    end

end

function flux = conserveFlux(state, G, bc)

f = G.cells.faces(:,1);
ncf = diff(G.cells.facePos);
fSgn = 1 - 2*(G.faces.neighbors(f,1) ~= rldecode((1:G.cells.num)', ncf,1));

[ii, jj] = blockDiagIndex(ones(G.cells.num,1), ncf);
r = -sparse(ii,jj,1)*(stateVEM.flux(f).*fSgn);

if norm(r) > 1e-15

    ncf = diff(G.cells.facePos);

    fn = bsxfun(@times,G.faces.normals(f,:), fSgn./G.faces.areas(f));
    delta = (fn*K)';
    [ii, jj] = blockDiagIndex(3*ones(numel(f),1), ones(numel(f),1));
    delta = sparse(ii,jj,delta(:));

    fn = fn';
    fn = sparse(ii,jj,fn(:));
    delta = fn'*delta;
    delta = delta*ones(size(delta,1),1);

    ii = f;
    jj = (1:numel(f))';
    P = sparse(ii, jj,1);
    omega = P*delta;

    for i = 1:G.faces.num
        d = delta(f == i);
        omega(i) = omega(i)/(numel(d)*prod(d));
    end

    B = zeros(G.cells.num, G.cells.num);

    for i = 1:G.cells.num
        fi = G.cells.faces(G.cells.facePos(i):G.cells.facePos(i+1)-1);
        for j = 1:G.cells.num
            fj = G.cells.faces(G.cells.facePos(j):G.cells.facePos(j+1)-1);
            ff = fi(ismember(fi, fj));
            ffSgn = 1-2*(G.faces.neighbors(ff,1) ~= i);
            B(i,j) = -sum(1./omega(ff).*G.faces.areas(ff).*ffSgn);
        end
    end

    beta = B\r;
    beta = rldecode(beta, ncf,1);
    I = sparse(f, 1:numel(f), 1);
    beta = I*beta.*G.faces.areas;

    flux = stateVEM.flux - 1./omega.*beta;

    [ii, jj] = blockDiagIndex(ones(G.cells.num,1), ncf);
    r = -sparse(ii,jj,1)*(stateVEM.flux(f).*fSgn);
    if norm(r)/norm(flux) > 1e-15;
        warning('Could not construct conservative flux field');
    end

else
    flux = state.flux;
end
    
end

%--------------------------------------------------------------------------

function totmob = dynamic_quantities(state, fluid)

   [mu, ~] = fluid.properties(state);
   s       = fluid.saturation(state);
   kr      = fluid.relperm(s, state);
   
   mob    = bsxfun(@rdivide, kr, mu);
   totmob = sum(mob, 2);

end

%--------------------------------------------------------------------------

function nk = polyDim(k, dim)

    if k < 0
        nk = 0;
    else
        nk = nchoosek(k+dim,k);
    end
end
