function mpfastruct = computeNeumannMultiPointTrans(G, rock, varargin)
% Compute multi-point transmissibilities for MPFA for Neumann boundary condition.
%
% SYNOPSIS:
%   function mpfastruct = computeNeumannMultiPointTrans(G, rock, varargin)
%
% DESCRIPTION:
%
% PARAMETERS:
%   G        - Grid
%   rock     - Rock data structure (see description in `PermTensor`)
%   varargin - see below
%
% KEYWORD ARGUMENTS:
%   verbose       - true if verbose
%   blocksize     - size of the blocks (last block will have different size to adjust to grid)
%   ip_compmethod - Option sent to blockLocalFluxMimeticAssembly
%   eta           - Option sent to blockLocalFluxMimeticAssembly
%   invertBlocks  - Method by which to invert a sequence of small matrices that
%                   arise in the discretisation.  String.  Must be one of
%                      - MATLAB -- Use an function implemented purely in MATLAB
%                                  (the default).
%    
%                      - MEX    -- Use two C-accelerated MEX functions to
%                                  extract and invert, respectively, the blocks
%                                  along the diagonal of a sparse matrix.  This
%                                  method is often faster by a significant
%                                  margin, but relies on being able to build
%                                  the required MEX functions.
% RETURNS:
%
%   mpfastruct with fields:
%
%               'A'   : system matrix 
%               'F'   : flux operator 
%               'tbls': table structure
%
% EXAMPLE:
%
% SEE ALSO:
%
%
%{
Copyright 2009-2020 SINTEF Digital, Mathematics & Cybernetics.

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


   opt = struct('verbose'      , mrstVerbose, ...
                'blocksize'    , []         , ...
                'ip_compmethod', 'general'  , ...
                'invertBlocks' , 'matlab'   , ...
                'eta'          , 0);

   % possible options for ip_compmethod
   % 'general'       : general case, no special requirements on corner
   % 'nicecorner'    : case where at each corner, number of faces is equal to G.griddim
   % 'directinverse' : case where at each corner, number of faces is equal to
   %                   G.griddim AND eta is value such that N*R is diagonal (see Lipnikov paper)
   
   opt = merge_options(opt, varargin{:});
   switch opt.ip_compmethod
     case {'general', 'nicecorner'}
     case 'directinverse'
       isOk = false;
       if (G.griddim == 2) & (opt.eta == 1/3), isOk = true, end
       if (G.griddim == 3) & (opt.eta == 1/4), isOk = true, end
       if ~isOk
           error(['option values for ip_compmethod and eta are not ' ...
                  'compatible']);
       end
     otherwise
       error('value of option nodeompcase is not recognized.');
   end
   
   opt.invertBlocks = blockInverter(opt);
   blocksize = opt.blocksize;
   
   if opt.verbose
       fprintf('Computing inner product on sub-half-faces ...\n');
       t0 = tic();
   end
   
   opts = {'verbose'      , opt.verbose      , ...
           'ip_compmethod', opt.ip_compmethod, ...
           'eta'          , opt.eta};
   
   [B, tbls] = computeLocalFluxMimetic(G, rock, opts{:});
   
   if opt.verbose
       t0 = toc(t0);
       fprintf('Computing inner product on sub-half-faces done in %g sec\n', t0);
   end
   
   %% Assemble of the divergence operator, from facenode values to cell value.
   nf = G.faces.num;
   celltbl = tbls.celltbl;
   facetbl = tbls.facetbl;
   cellnodefacetbl = tbls.cellnodefacetbl;
   facenodetbl = tbls.facenodetbl;
   facenodetbl = facenodetbl.addLocInd('fnind');
   
   extfaces = (G.faces.neighbors(:, 1) == 0) | (G.faces.neighbors(:, 2) == 0);
   intfaces = find(~extfaces);
   clear intfacetbl
   intfacetbl.faces = intfaces;
   intfacetbl = IndexArray(intfacetbl);
   
   % setup table with only internal faces
   cellintnodefacetbl = crossIndexArray(cellnodefacetbl, intfacetbl, {'faces'});
   % add facenode indexing
   cellintnodefacetbl = crossIndexArray(cellintnodefacetbl, facenodetbl, ...
                                        {'faces', 'nodes'});
   
   fino = cellintnodefacetbl.get('faces');
   cino = cellintnodefacetbl.get('cells');
   sgn = 2*(cino == G.faces.neighbors(fino, 1)) - 1;
   
   map = TensorMap();
   map.fromTbl = cellintnodefacetbl;
   map.toTbl = cellnodefacetbl;
   map.mergefds = {'cells', 'nodes', 'faces', 'cnfind'};
   map = map.setup();
   
   sgn = map.eval(sgn);
   
   prod = TensorProd();
   prod.tbl1 = cellnodefacetbl;
   prod.tbl2 = facenodetbl;
   prod.tbl3 = celltbl;
   prod.reducefds = {'nodes', 'faces'};
   prod = prod.setup();
   
   div_T = SparseTensor();
   div_T = div_T.setFromTensorProd(sgn, prod);
   div = div_T.getMatrix();
   
   %% Invert matrix B
   % The facenode degrees of freedom, as specified by the facenodetbl table, are
   % ordered by nodes first (see implementation below). It means in particular
   % that the matrix B is, by construction, block diagonal.

   if opt.verbose
       fprintf('Computing inverse mixed innerproduct...\n');
       t0 = tic();   
   end

   % if we know - a priori - that matrix is symmetric, then we remove
   % symmetry loss that has been introduced in assembly.
   if strcmp(opt.ip_compmethod, 'directinverse')
       B = 0.5*(B + B');
   end
   
   nodes = facenodetbl.get('nodes');
   [~, sz] = rlencode(nodes); 
   iB = opt.invertBlocks(B, sz);
   % if we know - a priori - that matrix is symmetric, then we remove the loss of
   % symmetry that may have been introduced by the numerical inversion.
   if strcmp(opt.ip_compmethod, 'directinverse')
       iB = 0.5*(iB + iB');
   end
   
   if opt.verbose
       t0 = toc(t0);   
       fprintf('... computing inverse mixed innerproduct done in %g sec\n', t0);
   end
   
   %% Assemble the system matrix operaror: The degrees of freedom are the pressure
   % values at the cell center and at the external facenode.
   A = div*iB*div';
   
   %% Assemble the flux operator: From pressure values at the cell center, compute the fluxes at the faces
   F = iB*div';
   % mapping from facenode to face
   intfacenodetbl = crossIndexArray(facenodetbl, intfacetbl, {'faces'});
   
   map = TensorMap();
   map.fromTbl = facenodetbl;
   map.toTbl = facetbl;
   map.mergefds = {'faces'};
   map = map.setup();
   
   S_T = SparseTensor();
   S_T = S_T.setFromTensorMap(map);
   S = S_T.getMatrix();
   
   F = S*F;

   mpfastruct = struct('div' , div , ...
                       'F'   , F   , ...
                       'A'   , A   , ...
                       'tbls', tbls);
   
end

