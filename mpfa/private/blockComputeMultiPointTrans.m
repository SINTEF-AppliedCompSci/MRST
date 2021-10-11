function mpfastruct = blockComputeMultiPointTrans(G, rock, varargin)
% Block assembly of the multipoint transmissibilities for MPFA
%
% SYNOPSIS:
%   function mpfastruct = blockComputeMultiPointTrans(G, rock, varargin)
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
%
% RETURNS:
%
%   mpfastruct with fields:
%
%               'iB'  : inverse of scalar product (facenode degrees of freedom)
%               'div' : divergence operator (facenode values to cell values)
%               'Pext': projection operator on external facenode values
%               'F'   : flux operator (from cell and external facenode values to facenode values)
%               'A'   : system matrix (cell and external facenode degrees of freedom)
%               'tbls': table structure
% EXAMPLE:
%
% SEE ALSO:
% `computeMultiPointTrans`, `private/computeMultiPointTransTensorAssembly`
%
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

   cellfacetbl.cells = rldecode((1 : nc)', diff(G.cells.facePos));
   cellfacetbl.faces = G.cells.faces(:, 1);
   cellfacetbl = IndexArray(cellfacetbl);

   % sign of normal (outwards or inwards)
   N = G.faces.neighbors;
   fino = cellfacetbl.get('faces');
   cino = cellfacetbl.get('cells');
   facesgn = 2*(cino == N(fino, 1)) - 1;

   facenodetbl.faces = rldecode((1 : nf)', diff(G.faces.nodePos));
   facenodetbl.nodes = G.faces.nodes;
   facenodetbl = IndexArray(facenodetbl);

   cellfacenodetbl = crossIndexArray(cellfacetbl, facenodetbl, {'faces'});

   cellnodetbl = projIndexArray(cellfacenodetbl, {'cells', 'nodes'});

   extfaces = (G.faces.neighbors(:, 1) == 0) | (G.faces.neighbors(:, 2) == 0);
   extfacetbl.faces = find(extfaces);
   extfacetbl = IndexArray(extfacetbl);
   extfacenodetbl = crossIndexArray(facenodetbl, extfacetbl, {'faces'});
   extfacenodetbl = extfacenodetbl.addLocInd('extfnind');

   next = extfacenodetbl.num;

   A11 = sparse(nc  , nc);
   A12 = sparse(nc  , next);
   A21 = sparse(next, nc);
   A22 = sparse(next, next);
   F1  = sparse(nf  , nc);
   F2  = sparse(nf  , next);

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
                                                 opt.ip_compmethod );
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

       ind = sub2ind([lfn_num, lfn_num], ind1, ind2);
       iB = iBmat(ind);

       map = TensorMap();
       map.fromTbl = cellfacetbl;
       map.toTbl = loccellfacenodetbl;
       map.mergefds = {'cells', 'faces'};
       map = map.setup();

       sgn = map.eval(facesgn);
       div = sgn; % same as operator

       locface2nodetbl = locface2nodetbl.duplicateInd({'nodes', {'nodes1', ...
                           'nodes2'}});

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

       % Aggregate contribution in A11
       cells1 = diviBdivtbl.get('cells1');
       cells2 = diviBdivtbl.get('cells2');
       locA11 = sparse(cells1, cells2, diviBdiv, nc, nc);
       A11 = A11 + locA11;

       % Aggregate contribution in F1 (used to compute flux)
       
       faceiBdivtbl = projIndexArray(iBdivtbl, {'cells2', 'faces1'});
       map = TensorMap();
       map.fromTbl = iBdivtbl;
       map.toTbl = faceiBdivtbl;
       map.mergefds = {'cells2', 'faces1'};
       map = map.setup;

       locF1 = map.eval(iBdiv);

       faces = faceiBdivtbl.get('faces1');
       cells = faceiBdivtbl.get('cells2');
       locF1 = sparse(faces, cells, locF1, nf, nc);
       F1 = F1 + locF1;

       % Assemble parts corresponding to the boundary conditions (if they exist for
       % this block)
       locfacetbl = projIndexArray(locfacenodetbl, {'faces'});
       locfaces = locfacetbl.get('faces');
       locextfaces = (G.faces.neighbors(locfaces, 1) == 0) | (G.faces.neighbors(locfaces, 2) == 0);

       if any(locextfaces)
           % In the following, the mapping Pext acts as a projection on the external
           % face-node degrees of freedom, but we also have to multiply the
           % coefficients by the sign of the normal to get the correct flux
           % direction.
           
           clear locextfacetbl
           locextfacetbl.faces = locfaces(locextfaces);
           locextfacetbl = IndexArray(locextfacetbl);

           locextfacenodetbl = crossIndexArray(locfacenodetbl, locextfacetbl, {'faces'});
           % we fetch index 'extfnind' for locextfacenodetbl from extfacenodetbl
           locextfacenodetbl = crossIndexArray(locextfacenodetbl, extfacenodetbl, {'faces', 'nodes'});

           map = TensorMap();
           map.fromTbl = loccellfacenodetbl;
           map.toTbl = locextfacenodetbl;
           map.mergefds = {'faces', 'nodes'};
           map = map.setup();

           Pext = map.eval(sgn);

           % We assemble the part corresponding to A21 = -Pext*iB*div';

           prod = TensorProd();
           prod.tbl1 = iBdivtbl;
           prod.tbl2 = locextfacenodetbl;
           prod.replacefds2 = {{'faces', 'faces1'}, ...
                               {'nodes', 'nodes1'}, ...
                               {'extfnind', 'extfnind1'}};
           prod.mergefds = {'faces1', 'nodes1'};
           prod = prod.setup();

           PextiBdiv = prod.eval(iBdiv, Pext);
           PextiBdivtbl = prod.tbl3;

           extfnind = PextiBdivtbl.get('extfnind1');
           cells = PextiBdivtbl.get('cells2');

           locA21 = sparse(extfnind, cells, PextiBdiv, next, nc);

           A21 = A21 - locA21;

           % We assemble the part corresponding to A12 = -div*iB*Pext';
           
           % First we compute  iB*Pext' (denoted iBPext).
           prod = TensorProd();
           prod.tbl1 = locface2nodetbl;
           prod.tbl2 = locextfacenodetbl;
           prod.replacefds2 = {{'faces', 'faces2'}, ...
                               {'nodes', 'nodes2'}, ...
                               {'extfnind', 'extfnind2'}};
           prod.mergefds = {'faces2', 'nodes2'};
           prod = prod.setup();

           iBPext = prod.eval(iB, Pext);
           iBPexttbl = prod.tbl3;

           % We multiply div with iBPexttbl 
           prod = TensorProd();
           prod.tbl1 = loccellfacenodetbl;
           prod.tbl2 = iBPexttbl;
           prod.replacefds1 = {{'cells', 'cells1'}, ...
                               {'faces', 'faces1'}, ...
                               {'nodes', 'nodes1'}};
           prod.reducefds = {'faces1', 'nodes1'};
           prod = prod.setup();

           diviBPext = prod.eval(div, iBPext);
           diviBPexttbl = prod.tbl3;
           
           extfnind = diviBPexttbl.get('extfnind2');
           cells = diviBPexttbl.get('cells1');

           locA12 = sparse(cells, extfnind, diviBPext, nc, next);
           A12 = A12 - locA12;
             
           % Aggregate contribution in F2 (contribution of external face-node to the face fluxes)
           
           locface_1node_2face_2tbl = projIndexArray(iBPexttbl, {'faces1', 'faces2', 'nodes2', 'extfnind2'});
           
           map = TensorMap();
           map.fromTbl = iBPexttbl;
           map.toTbl = locface_1node_2face_2tbl;
           map.mergefds = {'faces1', 'faces2', 'nodes2', 'extfnind2'};
           map = map.setup();

           faceiBPext = map.eval(iBPext);
           
           extfnind = locface_1node_2face_2tbl.get('extfnind2');
           faces = locface_1node_2face_2tbl.get('faces1');
           
           locF2 = sparse(faces, extfnind, faceiBPext, nf, next);
           
           F2 = F2 - locF2;

           % We assemble the part corresponding to A22 = Pext*iB*Pext';
           
           % we compute Pext*iBPext (denoted PextiBPext)
           
           prod = TensorProd();
           prod.tbl1 = locextfacenodetbl;
           prod.tbl2 = iBPexttbl;
           prod.replacefds1 = {{'faces', 'faces1'}, ...
                               {'nodes', 'nodes1'}, ...
                               {'extfnind', 'extfnind1'}};
           prod.mergefds = {'faces1', 'nodes1'};
           prod = prod.setup();
           
           PextiBPext = prod.eval(Pext, iBPext);
           PextiBPexttbl = prod.tbl3;
           
           extfnind1 = PextiBPexttbl.get('extfnind1');
           extfnind2 = PextiBPexttbl.get('extfnind2');
           
           locA22 = sparse(extfnind1, extfnind2, PextiBPext, next, next);
           A22 = A22 + locA22;
           
       end
       
       if opt.verbose
           fprintf('%g seconds\n', toc);
       end

   end

   tbls = struct('facenodetbl'    , facenodetbl    , ...
                 'cellfacenodetbl', cellfacenodetbl, ...
                 'extfacenodetbl' , extfacenodetbl);

   A = [[A11, A12]; [A21, A22]];
   F = [F1, F2];

   % filename = 'new';
   % savefiledebug;
   
   mpfastruct = struct('A'   , A, ...
                       'F'   , F, ...
                       'tbls', tbls);
end
