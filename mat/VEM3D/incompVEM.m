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
N  = G.nodes.num + G.faces.num*(k-1) + G.cells.num*k*(k-1)/2;

if G.griddim == 2
    nE = 0;
    nk = (k+1)*(k+2)/2;
    NK = diff(G.cells.nodePos) + diff(G.cells.facePos)*(k-1) + k*(k-1)/2;
else
    nE = G.edges.num;
    nk = (k+1)*(k+2)*(k+3)/6;
    NK = diff(G.cells.nodePos) + diff(G.cells.edgePos)*(k-1) ...
       + diff(G.cells.facePos)*k*(k-1)/2 + k*(k^2-1)/6;
end
nker = sum(NK - nk);

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

%%  ASSEMBLE GLOBAL MATRIX AND COMPUTE RIGHT-HAND SIDE                   %%

[A, rhs] = glob(G, S, opt.src, k, N);

[A, rhs] = imposeBC(G, opt.bc, k, A, rhs);

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

function [A, rhs] = glob(G, S, src, k, N)
    
    ii = 1:numel(S.dofVec); jj = S.dofVec;
    P = sparse(ii,jj,1, numel(S.dofVec), N);
    A = P'*S.A*P;

    if ~isempty(src)
        if k == 1
            rhs = zeros(G.cells.num,1);
            rhs(src.cell) = src.rate.*G.cells.volumes(src.cell);
            rhs = rldecode(rhs, diff(G.cells.nodePos), 1);
            rhs = S.PiN1'*S.PiN1*rhs;
            rhs = P'*rhs;
        else
            rhs = zeros(N,1);
            rhs(src.cell + G.nodes.num + G.faces.num) = src.rate';
        end
        
    else
        rhs = zeros(N,1);
    end
end

%--------------------------------------------------------------------------

function [A, rhs] = imposeBC(G, bc, k, A, rhs)

bc = interpolateBC(G, bc);

NK = G.nodes.num + G.faces.num*(k-1) + G.cells.num*k*(k-1)/2;

edges   = bc.face(strcmp(bc.type,'flux'));
edgeLengths = G.faces.areas(edges);

nodeNum = mcolon(G.faces.nodePos(edges), G.faces.nodePos(edges +1)-1);
nodes   = G.faces.nodes(nodeNum);
nN = numel(nodes);
uNodes = unique(nodes);
nUN = numel(uNodes);
S = (repmat(nodes,1,nUN) == repmat(uNodes',nN,1))';

vals = bc.value(strcmp(bc.type,'flux'),:);

if k == 1
    vals = vals.*[edgeLengths/6, edgeLengths/6, edgeLengths/3];
    vals = bsxfun(@plus, vals(:,1:2), vals(:,3));
    vals = reshape(vals',[],1);
    vals = S*vals;
    dofVec = uNodes';
elseif k == 2
    vals = vals.*[edgeLengths/6, edgeLengths/6, edgeLengths*2/3];
    vals = [S*reshape(vals(:,1:2)',[],1); vals(:,3)];
    dofVec = [uNodes', edges' + G.nodes.num];
end
rhs(dofVec) = rhs(dofVec) + vals;

edges   = bc.face(strcmp(bc.type,'pressure'));
nodeNum = mcolon(G.faces.nodePos(edges), G.faces.nodePos(edges +1)-1);
nodes   = G.faces.nodes(nodeNum);
vals = bc.value(strcmp(bc.type,'pressure'),:);

if k == 1
    vals = reshape(vals(:,1:2)',[],1);
    dofVec = nodes';
elseif k == 2
    dofVec = [nodes', edges' + G.nodes.num];
    vals = [reshape(vals(:,1:2)',[],1); vals(:,3)];
end
rhs(dofVec) = vals;
I = spdiags(ones(NK,1),0,NK,NK);
A(dofVec,:) = I(dofVec,:);

end

%--------------------------------------------------------------------------

function bc = interpolateBC(G, bc)

    value = zeros(numel(bc.face), 3);
    
    isNeu = strcmp(bc.type, 'flux');
    value(isNeu,:) = .5*repmat(bc.value(isNeu), 1, 3);
    
    isDir = strcmp(bc.type, 'pressure');
    f = bc.face(isDir);
    fSign = (-ones(numel(f),1)).^(G.faces.neighbors(f,2) ~= 0); 
    
    ve = bc.value(isDir);
    
    nodes = G.faces.nodes(mcolon(G.faces.nodePos(f), ...
                                 G.faces.nodePos(f+1)-1));
    if size(nodes,1) == 1; nodes = nodes'; end
    nodes   = reshape(nodes,2,[])';
    nodes(fSign == -1,:) = nodes(fSign == -1, 2:-1:1);
    
    
    M = bsxfun(@eq, repmat(nodes(:,2)', size(nodes,1), 1), nodes(:,1));
    pn = M*(1:numel(f))';
    nn = M'*(1:numel(f))';
        
    d = sqrt(sum((G.faces.centroids(f( pn(pn~= 0)),:)-G.faces.centroids(f(pn~=0),:)).^2,2));
    a = (ve(pn(pn ~= 0))-ve(pn ~= 0))./d;
    fa = G.faces.areas(f);
    
    vn = zeros(numel(nodes)/2 + nnz(pn==0),1);
    
    isInt = false(numel(vn),1);
    isInt(1) = pn(1);
    j = 2;
    for i = 2:numel(pn)
        isInt(j) = any([pn(i),nn(i-1)]~= 0, 2);
        j = j+1;
        if pn(i) == 0 && nn(i-1) == 0
            j = j+1;
        end
    end
    
    vn(isInt) = a.*fa(pn~=0)/2 + ve(pn~=0);
    
    ii = find(pn==0) - (0:nnz(pn==0)-1)';
    jj = find(pn==0) + (0:nnz(pn==0)-1)';
    vn(jj) = a(ii).*fa(pn==0)/2 + ve(pn==0);
    
    ii = find(nn==0) - (1:nnz(nn==0))';
    jj = find(nn==0) + (1:nnz(nn==0))';
    vn(jj) = -a(ii).*fa(nn==0)/2 + ve(nn==0);
    
    vn = reshape(rldecode(vn, isInt+1, 1),2,[])';
    value(isDir,:) = [vn, ve];
    
    bc.value = value;
      
end