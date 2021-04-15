function mpfastruct = blockComputeNeumannMultiPointTrans(G, rock, varargin)
% Block assembly of the multipoint transmissibilities for MPFA for Neumann boundary conditions
%
% SYNOPSIS:
%   function mpfastruct = blockComputeNeumannMultiPointTrans(G, rock, varargin)
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
%  `computeMultiPointTrans`, `private/computeMultiPointTransTensorAssembly`

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

   opt = struct('verbose'      , mrstVerbose, ...
                'blocksize'    , []         , ...
                'ip_compmethod', 'general'  , ...
                'invertBlocks' , 'matlab'   , ...
                'eta'          , 0);

   opt = merge_options(opt, varargin{:});
   opt.invertBlocks = blockInverter(opt);
   blocksize = opt.blocksize;
   
   nn = G.nodes.num;
   nblocks = floor(nn/blocksize);
   blocksizes = repmat(blocksize, nblocks, 1);
   if nn > nblocks*blocksize
       blocksizes = [blocksizes; nn - nblocks*blocksize];
   end
   nblocks = numel(blocksizes);
   blockinds = cumsum([1; blocksizes]);

   nc  = G.cells.num;
   nf  = G.faces.num;
   dim = G.griddim;
   
   N = G.faces.neighbors;
   
   facetbl.faces = (1 : nf)';
   facetbl = IndexArray(facetbl);
   
   % setup table of cell face index
   cellfacetbl.cells = rldecode((1 : nc)', diff(G.cells.facePos));
   cellfacetbl.faces = G.cells.faces(:, 1);
   cellfacetbl = IndexArray(cellfacetbl);

   % sign of normal (outwards or inwards)
   N = G.faces.neighbors;
   fino = cellfacetbl.get('faces');
   cino = cellfacetbl.get('cells');
   facesgn = 2*(cino == N(fino, 1)) - 1;
   
   
   % setup table of face node index
   facenodetbl.faces = rldecode((1 : nf)', diff(G.faces.nodePos));
   facenodetbl.nodes = G.faces.nodes;
   facenodetbl = IndexArray(facenodetbl);

   cellfacenodetbl = crossIndexArray(cellfacetbl, facenodetbl, {'faces'});
   
   % setup table of cell node index
   cellnodetbl = projIndexArray(cellfacenodetbl, {'cells', 'nodes'});      
   
   extfaces = (N(:, 1) == 0) | (N(:, 2) == 0);
   intfaces = ~extfaces;
   
   A = sparse(nc, nc);
   F = sparse(nf, nc);

   clear globtbls
   globtbls.cellfacetbl     = cellfacetbl;
   globtbls.facenodetbl     = facenodetbl;
   globtbls.cellfacenodetbl = cellfacenodetbl;
   globtbls.cellnodetbl     = cellnodetbl;
   
   for iblock = 1 : nblocks

       if opt.verbose
           fprintf('Starting with block %d/%d ... ', iblock, nblocks);
           tic
       end
       nodes = [blockinds(iblock) : (blockinds(iblock + 1) - 1)]';
       [B, tbls] = blockLocalFluxMimeticAssembly(G, rock, globtbls, nodes, ...
                                                 'eta', opt.eta, 'ip_compmethod', ...
                                                 opt.ip_compmethod);
       if isempty(B)
           % handle case when the nodes do not belong to any faces.
           if opt.verbose
               fprintf('%g seconds\n', toc);
           end
           break
       end

       locfacenodetbl     = tbls.facenodetbl;
       locface2nodetbl    = tbls.face2nodetbl;
       loccellfacenodetbl = tbls.cellfacenodetbl;

       loccellfacenodetbl = loccellfacenodetbl.removeInd({'cnfind'});
       locfacenodetbl     = locfacenodetbl.removeInd({'fnind'});
       locface2nodetbl    = locface2nodetbl.removeInd({'fnind1', 'fnind2'});
       
       % Assembly of B
       prod = TensorProd();
       prod.tbl1 = locface2nodetbl;
       prod.tbl2 = locfacenodetbl;
       prod.tbl3 = locfacenodetbl;
       prod.replacefds1 = {{'faces1', 'faces'}};
       prod.replacefds2 = {{'faces', 'faces2'}};
       prod.mergefds = {'nodes'};
       prod.reducefds = {'faces2'};
       prod = prod.setup();

       [ind1, ind2] = prod.getDispatchInd();
       lfn_num = locfacenodetbl.num;
       Bmat = sparse(ind1, ind2, B, lfn_num, lfn_num);
       
       % if we know - a priori - that matrix is symmetric, then we remove
       % symmetry loss that has been introduced in assembly.
       if strcmp(opt.ip_compmethod, 'directinverse')
           Bmat = 0.5*(Bmat + Bmat');
       end
       nodes = locfacenodetbl.get('nodes');
       [~, sz] = rlencode(nodes);
       iBmat   = opt.invertBlocks(Bmat, sz);
       % if we know - a priori - that matrix is symmetric, then we remove the loss of
       % symmetry that may have been induced by the numerical inversion.
       if strcmp(opt.ip_compmethod, 'directinverse')
           iBmat = 0.5*(iBmat + iBmat');
       end
       
       ind = sub2ind([lfn_num, lfn_num], ind2, ind1);
       iB = iBmat(ind);
       
       % clean-up and prepare locface2nodetbl for further use in contraction operations
       locface2nodetbl = locface2nodetbl.duplicateInd({'nodes', {'nodes1', ...
                           'nodes2'}});
       
       % remove external faces from loccellfacenodetbl and locfacenodetbl
       locfaces = loccellfacenodetbl.get('faces');
       isintfaces = (intfaces(locfaces) > 0);
       loccellfacenodetbl.inds = loccellfacenodetbl.inds(isintfaces, :);
       
       locfaces = locfacenodetbl.get('faces');
       isintfaces = (intfaces(locfaces) > 0);
       locfacenodetbl.inds = locfacenodetbl.inds(isintfaces, :);

       if locfacenodetbl.num == 0
           % handle case when all faces are external
           if opt.verbose
               fprintf('%g seconds\n', toc);
           end
           break  
       end
       
       map = TensorMap();
       map.fromTbl = cellfacetbl;
       map.toTbl = loccellfacenodetbl;
       map.mergefds = {'cells', 'faces'};
       map = map.setup();
       
       div = map.eval(facesgn);

       prod = TensorProd();
       prod.tbl1 = locface2nodetbl;
       prod.tbl2 = loccellfacenodetbl;
       prod.replacefds2 = {{'cells', 'cells2'}, {'faces', 'faces2'}, {'nodes', ...
                           'nodes2'}};
       prod.reducefds = {'faces2', 'nodes2'};
       prod = prod.setup();
       
       iBdiv = prod.eval(iB, div);
       iBdivtbl = prod.tbl3;
       
       prod = TensorProd();
       prod.tbl1 = loccellfacenodetbl;
       prod.tbl2 = iBdivtbl;
       prod.replacefds1 = {{'cells', 'cells1'}, {'faces', 'faces1'}, {'nodes', ...
                           'nodes1'}};
       prod.reducefds = {'faces1', 'nodes1'};
       prod = prod.setup();
       
       diviBdiv    = prod.eval(div, iBdiv);
       diviBdivtbl = prod.tbl3;

       % Aggregate contribution in A
       
       cells1 = diviBdivtbl.get('cells1');
       cells2 = diviBdivtbl.get('cells2');
       locA = sparse(cells1, cells2, diviBdiv, nc, nc);
       A = A + locA;
       
       % Aggregate contribution in F
       locface_1cell_2tbl = projIndexArray(iBdivtbl, {'faces1', 'cells2'});
       % remove external faces
       locfaces = locface_1cell_2tbl.get('faces1');
       isintfaces = (intfaces(locfaces) > 0);
       locface_1cell_2tbl.inds = locface_1cell_2tbl.inds(isintfaces, :);
       
       map = TensorMap();
       map.fromTbl = iBdivtbl;
       map.toTbl = locface_1cell_2tbl;
       map.mergefds = {'faces1', 'cells2'};
       map = map.setup;
       
       locF = map.eval(iBdiv);
       
       faces1 = locface_1cell_2tbl.get('faces1');
       cells2 = locface_1cell_2tbl.get('cells2');
       locF = sparse(faces1, cells2, locF, nf, nc);
       F    = F + locF;
       
       if opt.verbose
           fprintf('%g seconds\n', toc);
       end
   end

   tbls = struct('facenodetbl'    , facenodetbl, ...
                 'cellfacenodetbl', cellfacenodetbl);
   
   mpfastruct = struct('A'   , A   , ...
                       'F'   , F   , ...
                       'tbls', tbls);
end

