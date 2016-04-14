function [sol, varargout] = VEM2D(G, f, k, bc, varargin)
%--------------------------------------------------------------------------
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
%       or, if a fluid is specified,
%
%           -\Delta p = \frac{\mu}{\rho} f,
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
%       fluid        - Single phase fluid struct constructed using
%                      initSingleFluid.
%       projectors   - Boolean. If true, matrix representations
%                      of \Pi^\nabla in the monomial basis \mathcal_k(K)
%                      will be added to grid structure G.
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
%       [1]     - Thesis title.
%-----------------------------------------------------------------ØSK-2016-

%{
   Copyright (C) 2016 Øystein Strengehagen Klemetsdal. See Copyright.txt
   for details.
%}

addpath('./');

%%  MERGE INPUT PARAMETRES                                               %%

nN = G.nodes.num;
nE = G.faces.num;
nK = G.cells.num;

nk   = (k+1)*(k+2)/2;
NK   = diff(G.cells.nodePos) + diff(G.cells.facePos)*(k-1) + k*(k-1)/2;
nker = sum(NK - nk);

opt = struct('sigma'       , 1         , ...
             'src'         , []        , ...
             'fluid'       , []        , ...
             'projectors'  , false     , ...
             'cellAverages', false     );
opt = merge_options(opt, varargin{:});

sigma        = opt.sigma;
src          = opt.src;
fluid        = opt.fluid;
projectors   = opt.projectors;
cellAverages = opt.cellAverages;

if isempty(fluid)
    mu = 1; rho = 1;
else
    [mu, rho] = fluid.properties();
end

%%  CHECK CORRECTNESS OF INPUT                                           %%

assert(G.griddim == 2, 'VEM2D is only supproted for 2D grids');

if ~isa(f, 'function_handle')
    assert(numel(f) == 1, ...
             'Source function f must either be scalar or function handle');
end

assert(k == 1 | k == 2, 'VEM only implemented for 1st and second order');

assert(any(numel(sigma) == [sum(nker),1]), ...
     'Number of elements in paramter matrix sigma must be 1 or sum(nker)');

assert(islogical(projectors), ' ''projectors'' must be boolean')

assert(islogical(cellAverages), ' ''cellAverages'' must be boolean')
        
if cellAverages
    projectors = true;
end

%%  COMPUTE STIFFNESS MATRIX, LOAD TERM AND PROJECTORS                   %%

[A,b,PNstarT] = VEM2D_glob(G, f, k, bc, sigma, projectors, src, mu, rho);

%%  SOLVE LINEAR SYSTEM                                                  %%

fprintf('Solving linear system ...\n')
tic;

U = A\b;

stop = toc;
fprintf('Done in %f seconds.\n\n', stop);

%%  MAKE SOLUTION STRUCT                                                 %%

nodeValues  = full( U( 1:nN)                            );
edgeValues  = full( U((1:nE*(k-1)) + nN)                );
cellMoments = full( U((1:nK*k*(k-1)/2) + nN + nE*(k-1)) );

sol = struct(...
             'nodeValues' , {nodeValues} , ...
             'edgeValues' , {edgeValues} , ...
             'cellMoments', {cellMoments}     );
if projectors
    G.cells.('PNstarT') = PNstarT;
    PNstarPos = [1, cumsum(diff(G.cells.nodePos') + ...
                           diff(G.cells.facePos')*(k-1) + k*(k-1)/2) + 1];
    G.cells.('PNstarPos') = PNstarPos;
    varargout(1) = {G};
end

if cellAverages && k == 1
    sol = calculateCellAverages(G, sol);
end
         
end