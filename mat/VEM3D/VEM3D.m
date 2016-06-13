function [sol, varargout] = VEM3D(G, f, bc, k, varargin)
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

addpath('../');

%%  MERGE INPUT PARAMETRES                                               %%

nN = G.nodes.num;
nE = G.edges.num;
nF = G.faces.num;
nK = G.cells.num;

nk   = (k+1)*(k+2)*(k+3)/6;
NK   = diff(G.cells.nodePos) + diff(G.cells.edgePos)*(k-1) ...
       + diff(G.cells.facePos)*k*(k-1)/2 + k*(k^2-1)/6;
nker = sum(NK - nk);

opt = struct('sigma'          , 1     , ...
             'src'            , []    , ...
             'faceProjectors' , false , ...
             'cellProjectors' , false , ...
             'faceAverages'   , false , ...
             'cellAverages'   , false     );
opt = merge_options(opt, varargin{:});

sigma          = opt.sigma;
src            = opt.src;
faceProjectors = opt.faceProjectors;
cellProjectors = opt.cellProjectors;
faceAverages   = opt.faceAverages;
cellAverages   = opt.cellAverages;

%%  CHECK CORRECTNESS OF INPUT                                           %%

assert(G.griddim == 3, 'VEM3D is only supported for 3D grids');

if ~isa(f, 'function_handle')
    assert(numel(f) == 1, ...
             'Source function f must either be scalar or function handle');
end

assert(k == 1 | k == 2, 'VEM only implemented for 1st and 2nd order');

assert(any(numel(sigma) == [sum(nker),1]), ...
     'Number of elements in paramter matrix sigma must be 1 or sum(nker)');

assert(islogical(cellProjectors), ' ''cellProjectors'' must be boolean')

assert(islogical(faceAverages), ' ''faceAverages'' must be boolean')

assert(islogical(cellAverages), ' ''cellAverages'' must be boolean')
        
if cellAverages
    cellProjectors = true;
end

%%  COMPUTE STIFFNESS MATRIX, LOAD TERM AND PROJECTORS                   %%

[A,b,G] ...
            = VEM3D_glob(G, f, bc, k, sigma, cellProjectors, src);

%%  SOLVE LINEAR SYSTEM                                                  %%

fprintf('Solving linear system ...\n')
tic;

U = A\b;

stop = toc;
fprintf('Done in %f seconds.\n\n', stop);

%%  MAKE SOLUTION STRUCT                                                 %%

nodeValues  = full( U( 1:nN)                                             );
edgeValues  = full( U((1:nE*(k-1))       + nN)                           );
faceMoments = full( U((1:nF*k*(k-1)/2)   + nN + nE*(k-1))                );
cellMoments = full( U((1:nK*k*(k^2-1)/6) + nN + nE*(k-1) + nF*k*(k-1)/2) );

sol = struct(...
             'nodeValues' , {nodeValues} , ...
             'edgeValues' , {edgeValues} , ...
             'faceMoments', {faceMoments}, ...
             'cellMoments', {cellMoments}     );

if any([faceProjectors, cellProjectors])
    varargout(1) = {G};
end

if faceAverages && k == 1
    sol = calculateFaceAverages(G, sol);
end

if cellAverages && k == 1
    sol = calculateCellAverages(G, sol);
end
         
end