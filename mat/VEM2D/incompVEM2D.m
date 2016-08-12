function [state, varargout] = incompVEM2D(state, G, rock, fluid, k, varargin)
%   stateves the 2D Poisson equation using a kth order virtual element
%   method.
%
%   SYNOPSIS:
%       [state, varargout] = VEM2D(G, f, k, bc, varargin)
%
%   DESCRIPTION:
%       stateves the Poisson equation
%
%           -\Delta u = f
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
%       cartGridQ    - If ture, and G is a cartesian grid, the matrices Q
%                      in the stability term are set to the degrees of
%                      freedom of \sqrt(9/(4 h_x h_y) xy. See [1].
%       src          - Source term struct constructed using addSource.
%       projectors   - Boolean. If true, matrix representations
%                      of \Pi^\nabla in the monomial basis \mathcal_k(K)
%                      will be added to grid structure G.
%       cellAverages - Boolean. If true, exact cell averages of
%                      approximated stateution will be calculated
%                      for 1st order VEM. Useful for countour plots.
%
%   RETURNS:
%       state          - stateution struct. Contans the fileds
%                           * nodeValues, values at the nodes.
%                           * edgeValues, values at the edge
%                             midpoints. Empty for k = 1.
%                           * cellMoments, the first moment (avearge) over
%                             each cell. Empty for k = 1 unless
%                             cellAverages = true.
%
%   OPTIONAL RETURN VALUE:
%       G            - If projectors = true or cellAverages = true, grid
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
%       state  = VEM2D(G,f,bc,2);
%
%   REFERENCES:
%       [1]     - The virtual element method as a common framework for
%                 finite element and finite difference methods - Numerical
%                 and theoretical analysis.
%-----------------------------------------------------------------ØSK-2016-

%{
   Copyright (C) 2016 Øystein Strengehagen Klemetsdal. See Copyright.txt
   for details.
%}

%%  MERGE INPUT PARAMETRES                                               %%

nN = G.nodes.num;
nE = G.faces.num;
nK = G.cells.num;

nk   = (k+1)*(k+2)/2;
NK   = diff(G.cells.nodePos) + diff(G.cells.facePos)*(k-1) + k*(k-1)/2;
nker = sum(NK - nk);

opt = struct('bc'          , []         , ...
             'src'         , []         , ...
             'srcFunc'     , []         , ...
             'sigma'       , 1          , ...
             'cartGridQ'   , false      , ...
             'projectors'  , false      , ...
             'cellPressure', false      , ...
             'linSolve'    , @mldivide  , ...
             'matrixOutput', false           );
opt = merge_options(opt, varargin{:});

%%  CHECK CORRECTNESS OF INPUT                                           %%

assert(G.griddim == 2, 'VEM2D is only supproted for 2D grids')

if isempty(opt.srcFunc)
    opt.srcFunc = 0;
end

if ~isa(opt.srcFunc, 'function_handle')
    assert(numel(opt.srcFunc) == 1 || 0, ...
    'Source function ''srcFunc'' must either be scalar or function handle')
end

assert(k == 1 | k == 2, 'VEM only implemented for 1st and 2nd order')

assert(any(numel(opt.sigma) == [sum(nker),1]), ...
    'Number of elements in parameter matrix sigma must be 1 or sum(nker)')

assert(islogical(opt.projectors), ' ''projectors'' must be boolean')

assert(islogical(opt.cellPressure), ' ''cellPressure'' must be boolean')
        
if opt.cellPressure
    opt.projectors = true;
end

if isempty(opt.bc)
    opt.bc = VEM2D_addBC(opt.bc, G, boundaryFaces(G), 'flux', 0);
end

%%  COMPUTE STIFFNESS MATRIX, LOAD TERM AND PROJECTORS                   %%

[A,b,PNstarT] = VEM2D_glob(G, rock, fluid, k, opt.bc, opt.src, ...
                    opt.srcFunc, opt.sigma, opt.cartGridQ, opt.projectors);

%%  SOLVE LINEAR SYSTEM                                                  %%

fprintf('Solving linear system ...\n')
tic;

U = opt.linSolve(A,b);

stop = toc;
fprintf('Done in %f seconds.\n\n', stop);

if opt.matrixOutput
    state.A   = A;
    state.rhs = b;
end

%%  MAKE SOLUTION STRUCT                                                 %%

state.nodePressure = full( U( 1:nN)                            );
state.facePressure = full( U((1:nE*(k-1)) + nN)                );
state.pressure     = full( U((1:nK*k*(k-1)/2) + nN + nE*(k-1)) );

if opt.projectors
    G.cells.('PNstarT') = PNstarT;
    PNstarPos = [1, cumsum(diff(G.cells.nodePos') + ...
                           diff(G.cells.facePos')*(k-1) + k*(k-1)/2) + 1];
    G.cells.('PNstarPos') = PNstarPos;
    varargout(1) = {G};
end

if opt.cellPressure && k == 1
    state.pressure = calculatePressure(G, state);
end
         
end