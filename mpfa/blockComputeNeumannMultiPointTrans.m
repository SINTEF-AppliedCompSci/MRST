function mpfastruct = blockComputeNeumannMultiPointTrans(G, rock, varargin)
% Compute multi-point transmissibilities by block around nodes.


%{
Copyright 2009-2018 SINTEF Digital, Mathematics & Cybernetics.

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
   faces = (1 : nf)';
   extfaces = (N(:, 1) == 0) | (N(:, 2) == 0);
   intfaces = find(~extfaces);
   intfacetbl.faces = intfaces;
   intfacetbl.num = numel(intfacetbl.faces);
   
   cellfacetbl.cells = rldecode((1 : nc)', diff(G.cells.facePos));
   cellfacetbl.faces = G.cells.faces(:, 1);
   cellfacetbl.num   = numel(cellfacetbl.cells);

   [~, cellfacetbl] = setupTableMapping(cellfacetbl, intfacetbl, {'faces'});
   
   facenodetbl.faces = rldecode((1 : nf)', diff(G.faces.nodePos));
   facenodetbl.nodes = G.faces.nodes;
   facenodetbl.num = numel(facenodetbl.faces);

   [~, facenodetbl] = setupTableMapping(facenodetbl, intfacetbl, {'faces'});
   
   [~, cellfacenodetbl] = setupTableMapping(cellfacetbl, facenodetbl, ...
                                                         {'faces'});
   A = sparse(nc, nc);
   F = sparse(nf, nc);
   
   for iblock = 1 : nblocks

       if opt.verbose
           fprintf('Starting with block %d/%d ... ', iblock, nblocks);
           tic
       end
       nodes = [blockinds(iblock) : (blockinds(iblock + 1) - 1)]';
       [B, tbls] = blockLocalFluxMimeticAssembly(G, rock, nodes, 'eta', opt.eta, ...
                                                 'ip_compmethod', opt.ip_compmethod);
       
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

       loccellfacenodetbl = rmfield(loccellfacenodetbl, 'cnfind');
       locfacenodetbl = rmfield(locfacenodetbl, 'fnind');
       
       % Assembly of B
       Bmat = sparse(locface2nodetbl.fnind1, locface2nodetbl.fnind2, B, ...
                     locfacenodetbl.num, locfacenodetbl.num);

       [~, intlocfacenodetbl] = setupTableMapping(locfacenodetbl, intfacetbl, ...
                                                                {'faces'});
       if intlocfacenodetbl.num == 0
           if opt.verbose
               fprintf('%g seconds\n', toc);
           end
           break
       end
       
       % we order by node face to get block diagonal structure of the
       % resulting matrix
       intlocfacenodetbl = sortTable(intlocfacenodetbl, {'nodes', 'faces'});
       

       
       map = setupTableMapping(intlocfacenodetbl, locfacenodetbl, {'faces', ...
                           'nodes'});
       Bmat = map'*Bmat*map;
       [~, sz] = rlencode(intlocfacenodetbl.nodes);
       iBmat   = opt.invertBlocks(Bmat, sz);
       
       % remove from local table the external faces
       locfacenodetbl = intlocfacenodetbl;
       [~, loccellfacenodetbl] = setupTableMapping(loccellfacenodetbl, intfacetbl, ...
                                                                 {'faces'});
       locfacenodetbl = addLocInd(locfacenodetbl, 'fnind');
       % generate locface2nodetbl
       [~, locface2nodetbl] = setupTableMapping(locfacenodetbl, locfacenodetbl, ...
                                                              {'nodes'}, ...
                                                              'crossextend', ...
                                                              {{'faces', ...
                           {'faces1', 'faces2'}}, {'fnind', {'fnind1', ...
                           'fnind2'}}});
       
       clear locmattbl
       locmattbl.fnind = (1 : locfacenodetbl.num)';
       locmattbl.num   = locfacenodetbl.num;
       [~, locmattbl] = setupTableMapping(locmattbl, locmattbl, [], 'crossextend', ...
                                                     {{'fnind', {'fnind1', ...
                           'fnind2'}}});
       locmattbl = sortTable(locmattbl, {'fnind2', 'fnind1'});
       
       iB = iBmat(:);
       map = setupTableMapping(locmattbl, locface2nodetbl, {'fnind1', 'fnind2'});
       iB = map*iB;
       clear map;
       
       % clean-up and prepare locface2nodetbl for further use in contraction operations
       locface2nodetbl = duplicatefield(locface2nodetbl, {'nodes', {'nodes1', ...
                           'nodes2'}});
       locface2nodetbl = rmfield(locface2nodetbl, 'fnind1');
       locface2nodetbl = rmfield(locface2nodetbl, 'fnind2');
       
       div        = zeros(loccellfacenodetbl.num, 1);
       locfacetbl = projTable(locfacenodetbl, {'faces'});
       locfaces   = locfacetbl.faces;

       intn = (N(locfaces, 1) > 0);
       if any(intn)
           clear locposcellfacetbl
           locposcellfacetbl.cells = N(locfaces(intn), 1);
           locposcellfacetbl.faces = locfaces(intn);
           locposcellfacetbl.num   = numel(locposcellfacetbl.cells);
           mappos = setupTableMapping(locposcellfacetbl, loccellfacenodetbl, {'cells', 'faces'});
           div = div + mappos*ones(locposcellfacetbl.num, 1);
           clear mappos
       end

       intn = (N(locfaces, 2) > 0);
       if any(intn)
           clear locnegcellfacetbl
           locnegcellfacetbl.cells = N(locfaces(intn), 2);
           locnegcellfacetbl.faces = locfaces(intn);
           locnegcellfacetbl.num   = numel(locnegcellfacetbl.cells);
           mapneg = setupTableMapping(locnegcellfacetbl, loccellfacenodetbl, {'cells', 'faces'});
           div = div - mapneg*ones(locnegcellfacetbl.num, 1);
           clear mapneg
       end

       % Table loccell_1facenode_1tbl for div mapping from facenode to cell
       loccell_1facenode_1tbl = replacefield(loccellfacenodetbl, {{'faces', ...
                           'faces1'}, {'cells', 'cells1'}, {'nodes', ...
                           'nodes1'}});

       % Table loccell_2facenode_2tbl for div' mapping from cell to facenode
       loccell_2facenode_2tbl = replacefield(loccellfacenodetbl, {{'faces', ...
                           'faces2'}, {'cells', 'cells2'}, {'nodes', ...
                           'nodes2'}});
       
       [iBdiv, locfacenode_1cell_2tbl] = contractTable({iB, locface2nodetbl}, ...
                                                       {div, loccell_2facenode_2tbl}, ...
                                                       {{'nodes1', 'faces1'}, ...
                           {'cells2'}, {'faces2', 'nodes2'}});


       [diviBdiv, loccell_1cell_2tbl] = contractTable({div, loccell_1facenode_1tbl}, ...
                                                      {iBdiv, ...
                           locfacenode_1cell_2tbl}, {{'cells1'}, {'cells2'}, ...
                           {'faces1', 'nodes1'}});

       % Aggregate contribution in A
       tbl = loccell_1cell_2tbl; %alias
       locA = sparse(tbl.cells1, tbl.cells2, diviBdiv, nc, nc);
       A = A + locA;
       
       locface_1cell_2tbl = projTable(locfacenode_1cell_2tbl, {'faces1', ...
                           'cells2'});
       map = setupTableMapping(locfacenode_1cell_2tbl, locface_1cell_2tbl, ...
                                             {'faces1', 'cells2'});
       locF = map*iBdiv;
       tbl  = locface_1cell_2tbl; %alias
       locF = sparse(tbl.faces1, tbl.cells2, locF, nf, nc);
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

