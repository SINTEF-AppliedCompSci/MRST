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

if isempty(opt.bc)
    if G.griddim == 2
        opt.bc = VEM2D_addBC(opt.bc, G, boundaryFaces(G), 'flux', 0);
    else
        opt.bc = VEM3D_addBC(opt.bc, boundaryFaces(G), 'flux', 0);
    end
end

if opt.cellPressure
    opt.cellProjectors = true;
end


[mu, rho] = fluid.properties();
%%  ASSEMBLE GLOBAL MATRIX AND COMPUTE RIGHT-HAND SIDE                   %%

[A, rhs] = glob(G, S, opt.src, k, N, mu);

[A, rhs] = imposeBC(G, opt.bc, k, N, mu, A, rhs);

%%  SOLVE LINEAR SYSTEM                                                  %%

U = opt.linSolve(A, rhs);

%%  UPDATE STATE                                                         %%

state.nodePressure = ...
              full( U( 1:nN)                                             );
state.edgePressure = ...
              full( U((1:nE*(k-1))       + nN)                           );
state.facePressure = ...
              full( U((1:nF*k*(k-1)/2)   + nN + nE*(k-1))                );
state.cellPressure = ...
              full( U((1:nP*k*(k^2-1)/6) + nN + nE*(k-1) + nF*k*(k-1)/2) );

if any([opt.faceProjectors, opt.cellProjectors])
    varargout(1) = {G};
end

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

end

%--------------------------------------------------------------------------

function [A, rhs] = glob(G, S, src, k, N, mu)
    
        P = sparse(1:numel(S.dofVec), S.dofVec, 1, numel(S.dofVec), N);
    A = P'*S.A*P;

    if ~isempty(src)
        if k == 1
            rhs = zeros(G.cells.num,1);
            rhs(src.cell) = src.rate.*G.cells.volumes(src.cell);
            rhs = rldecode(rhs, diff(G.cells.nodePos), 1);
            rhs = S.PiN1'*S.PiN1*rhs;
            rhs = mu*P'*rhs;
        else
            rhs = zeros(N,1);
            rhs(src.cell + G.nodes.num + G.faces.num) = mu*src.rate';
        end
        
    else
        rhs = zeros(N,1);
    end
end

%--------------------------------------------------------------------------

function [A, rhs] = imposeBC(G, bc, k, N, mu, A, rhs)

if G.griddim == 2

    if ~isfield(bc, 'func') 
        bc = mrst2vem(bc, G);
    else
        bc.value = zeros(numel(bc.face), 3);
        for i = 1:numel(bc.face)
            g = bc.func{i};
            n = G.faces.nodes(G.faces.nodePos(bc.face(i)):G.faces.nodePos(bc.face(i)+1)-1);
            bc.value(i,:) = [g(G.nodes.coords(n,:))', ...
                             g(G.faces.centroids(bc.face(i),:))];
        end
    end

    f   = bc.face(strcmp(bc.type,'flux'));
    fa = G.faces.areas(f);

    n   = G.faces.nodes(mcolon(G.faces.nodePos(f), G.faces.nodePos(f +1)-1));
    nn = numel(n);
    un = unique(n);
    nun = numel(un);
    S = (repmat(n,1,nun) == repmat(un',nn,1))';

    v = mu*bc.value(strcmp(bc.type,'flux'),:);

    if k == 1
        v = v.*[fa/6, fa/6, fa/3];
        v = bsxfun(@plus, v(:,1:2), v(:,3));
        v = reshape(v',[],1);
        v = S*v;
        dofVec = un';
    elseif k == 2
        v = v.*[fa/6, fa/6, fa*2/3];
        v = [S*reshape(v(:,1:2)',[],1); v(:,3)];
        dofVec = [un', f' + G.nodes.num];
    end
    rhs(dofVec) = rhs(dofVec) + v;

    f   = bc.face(strcmp(bc.type,'pressure'));
    n   = G.faces.nodes(mcolon(G.faces.nodePos(f), G.faces.nodePos(f +1)-1));
    v = bc.value(strcmp(bc.type,'pressure'),:);

    if k == 1
        v = reshape(v(:,1:2)',[],1);
        dofVec = n';
    elseif k == 2
        dofVec = [n', f' + G.nodes.num];
        v = [reshape(v(:,1:2)',[],1); v(:,3)];
    end
    rhs(dofVec) = v;
    I = spdiags(ones(N,1),0,N,N);
    A(dofVec,:) = I(dofVec,:);
else

    isDir = strcmp(bc.type, 'pressure');
    f = bc.face(isDir);


    for i = 1:numel(f)
        n = G.faces.nodes(G.faces.nodePos(f(i)):G.faces.nodePos(f(i)+1)-1);
        e = G.faces.edges(G.faces.edgePos(f(i)):G.faces.edgePos(f(i)+1)-1);
        g = bc.func{i};
        rhs([n; e + G.nodes.num; f + G.nodes.num + G.edges.num]) ...
            = [g(G.nodes.coords(n,:)); g(G.edges.centroids(e,:)); ...
               polygonInt3D(G, f, g, k+1)./G.faces.areas(f)];
    end

    n = G.faces.nodes(mcolon(G.faces.nodePos(f), G.faces.nodePos(f+1)-1));
    e = G.faces.edges(mcolon(G.faces.edgePos(f), G.faces.edgePos(f+1)-1));
    dofVec = [n; e + G.nodes.num; f + G.nodes.num + G.edges.num];
    
    I = spdiags(ones(N,1),0,N,N);
    A(dofVec,:) = I(dofVec, :);
    
end

end

%--------------------------------------------------------------------------

function bc = mrst2vem(bc, G)
    
    isNeu = strcmp(bc.type, 'flux');
    
    bc.value(isNeu) ...
        = bc.value(isNeu)./G.faces.areas(bc.face(isNeu));

    bc.value = repmat(bc.value, 1,3);
    
    isDir = strcmp(bc.type, 'pressure');
    
    gD = scatteredInterpolant(G.faces.centroids(bc.face(isDir),:), ...
                                                        bc.value(isDir,1));
    n = G.faces.nodes(mcolon(G.faces.nodePos(bc.face(isDir)), ...
                             G.faces.nodePos(bc.face(isDir)+1)-1));
    vn = gD(G.nodes.coords(n,:));
    
    bc.value(isDir,1:2) = reshape(vn,2,[])';
    
end

%--------------------------------------------------------------------------

function nk = polyDim(k, dim)
    nk = nchoosek(k+dim,k);
end
