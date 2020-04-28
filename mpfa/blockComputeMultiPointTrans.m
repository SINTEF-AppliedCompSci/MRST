function mpfastruct = blockComputeMultiPointTrans(G, rock, varargin)
%Compute multi-point transmissibilities.
%
% SYNOPSIS:
%   T = computeMultiPointTrans2(G, rock)
%
% REQUIRED PARAMETERS:
%   G       - Grid data structure as described by grid_structure.
%
%   rock    - Rock data structure with valid field 'perm'.  The
%             permeability is assumed to be in measured in units of
%             metres squared (m^2).  Use function 'darcy' to convert from
%             (milli)darcies to m^2, e.g.,
%
%                 perm = convertFrom(perm, milli*darcy)
%
%             if the permeability is provided in units of millidarcies.
%
%             The field rock.perm may have ONE column for a scalar
%             permeability in each cell, TWO/THREE columns for a diagonal
%             permeability in each cell (in 2/3 D) and THREE/SIX columns
%             for a symmetric full tensor permeability.  In the latter
%             case, each cell gets the permeability tensor
%
%                 K_i = [ k1  k2 ]      in two space dimensions
%                       [ k2  k3 ]
%
%                 K_i = [ k1  k2  k3 ]  in three space dimensions
%                       [ k2  k4  k5 ]
%                       [ k3  k5  k6 ]
%
% OPTIONAL PARAMETERS
%   verbose   - Whether or not to emit informational messages throughout the
%               computational process.  Default value depending on the
%               settings of function 'mrstVerbose'.
%
%   invertBlocks -
%               Method by which to invert a sequence of small matrices that
%               arise in the discretisation.  String.  Must be one of
%                  - MATLAB -- Use an function implemented purely in MATLAB
%                              (the default).
%
%                  - MEX    -- Use two C-accelerated MEX functions to
%                              extract and invert, respectively, the blocks
%                              along the diagonal of a sparse matrix.  This
%                              method is often faster by a significant
%                              margin, but relies on being able to build
%                              the required MEX functions.
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

% COMMENTS:
%   PLEASE NOTE: Face normals have length equal to face areas.
%
% SEE ALSO:
%   `incompMPFAbc`, `mrstVerbose`.

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
   
   facetbl.faces = (1 : nf)';
   facetbl.num = nf;
   
   cellfacetbl.cells = rldecode((1 : nc)', diff(G.cells.facePos));
   cellfacetbl.faces = G.cells.faces(:, 1);
   cellfacetbl.num   = numel(cellfacetbl.cells);

   facenodetbl.faces = rldecode((1 : nf)', diff(G.faces.nodePos));
   facenodetbl.nodes = G.faces.nodes;
   facenodetbl.num = numel(facenodetbl.faces);

   [~, cellfacenodetbl] = setupTableMapping(cellfacetbl, facenodetbl, ...
                                                         {'faces'});
   
   cellnodetbl = projTable(cellfacenodetbl, {'cells', 'nodes'});   
  
   extfaces = (G.faces.neighbors(:, 1) == 0) | (G.faces.neighbors(:, 2) == 0);
   extfacetbl.faces = find(extfaces);
   extfacetbl.num   = numel(extfacetbl.faces);
   [~, extfacenodetbl] = setupTableMapping(facenodetbl, extfacetbl, {'faces'});
   extfacenodetbl = addLocInd(extfacenodetbl, 'extfnind');
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

       loccellfacenodetbl = rmfield(loccellfacenodetbl, 'cnfind');
       locfacenodetbl = rmfield(locfacenodetbl, 'fnind');
       
       % Assembly of B
       Bmat = sparse(locface2nodetbl.fnind1, locface2nodetbl.fnind2, B, ...
                     locfacenodetbl.num, locfacenodetbl.num);
       % if we know - a priori - that matrix is symmetric, then we remove
       % symmetry loss that has been introduced in assembly.
       if strcmp(opt.ip_compmethod, 'directinverse')
           Bmat = 0.5*(Bmat + Bmat');
       end
       [~, sz] = rlencode(locfacenodetbl.nodes);
       iBmat   = opt.invertBlocks(Bmat, sz);
       % if we know - a priori - that matrix is symmetric, then we remove the loss of
       % symmetry that may have been induced by the numerical inversion.
       if strcmp(opt.ip_compmethod, 'directinverse')
           iBmat = 0.5*(iBmat + iBmat');
       end
       [fnind1, fnind2, iB] = find(iBmat);
       clear locmattbl
       locmattbl.fnind1 = fnind1;
       locmattbl.fnind2 = fnind2;
       locmattbl.num = numel(locmattbl.fnind1);
       map = setupTableMapping(locmattbl, locface2nodetbl, {'fnind1', 'fnind2'});
       iB = map*iB;
       clear map;

       sgn        = zeros(loccellfacenodetbl.num, 1);
       locfacetbl = projTable(locfacenodetbl, {'faces'});
       locfaces   = locfacetbl.faces;

       intn = (N(locfaces, 1) > 0);
       if any(intn)
           clear locposcellfacetbl
           locposcellfacetbl.cells = N(locfaces(intn), 1);
           locposcellfacetbl.faces = locfaces(intn);
           locposcellfacetbl.num   = numel(locposcellfacetbl.cells);
           mappos = setupTableMapping(locposcellfacetbl, loccellfacenodetbl, {'cells', 'faces'});
           sgn = sgn + mappos*ones(locposcellfacetbl.num, 1);
           clear mappos
       end

       intn = (N(locfaces, 2) > 0);
       if any(intn)
           clear locnegcellfacetbl
           locnegcellfacetbl.cells = N(locfaces(intn), 2);
           locnegcellfacetbl.faces = locfaces(intn);
           locnegcellfacetbl.num   = numel(locnegcellfacetbl.cells);
           mapneg = setupTableMapping(locnegcellfacetbl, loccellfacenodetbl, {'cells', 'faces'});
           sgn = sgn - mapneg*ones(locnegcellfacetbl.num, 1);
           clear mapneg
       end
       div = sgn; % the two are equivalent as "operators"
       
       locface2nodetbl = locface2nodetbl.duplicateInd({'nodes', {'nodes1', ...
                           'nodes2'}});
       locface2nodetbl = rmfield(locface2nodetbl, 'fnind1');
       locface2nodetbl = rmfield(locface2nodetbl, 'fnind2');
       
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
       % Aggregate contribution in A11
       tbl = loccell_1cell_2tbl; %alias
       locA11 = sparse(tbl.cells1, tbl.cells2, diviBdiv, nc, nc);
       A11 = A11 + locA11;
       % Aggregate contribution in F1
       locface_1cell_2tbl = projTable(locfacenode_1cell_2tbl, {'faces1', ...
                           'cells2'});
       map = setupTableMapping(locfacenode_1cell_2tbl, locface_1cell_2tbl, ...
                                             {'faces1', 'cells2'});
       locF1 = map*iBdiv;
       tbl = locface_1cell_2tbl; %alias
       locF1 = sparse(tbl.faces1, tbl.cells2, locF1, nf, nc);
       F1 = F1 + locF1;
       
       % Assemble parts corresponding to the boundary conditions
       
       locfacetbl = projTable(locfacenodetbl, {'faces'});
       locfaces = locfacetbl.faces;
       locextfaces = (G.faces.neighbors(locfaces, 1) == 0) | (G.faces.neighbors(locfaces, 2) == 0);
       % if there are any boundary in this block partition
       if any(locextfaces)
           clear locfacexttbl
           locextfacetbl.faces = locfaces(locextfaces);
           locextfacetbl.num   = numel(locextfacetbl.faces);
           [~, locextfacenodetbl] = setupTableMapping(locfacenodetbl, locextfacetbl, ...
                                                                    {'faces'});
           op = setupTableMapping(loccellfacenodetbl, locextfacenodetbl, {'faces', 'nodes'});
           extsgn = op*sgn;
           
           % Table locextfacenode_1facenode_1tbl for Pext
           locextfacenode_1facenode_1tbl = locextfacenodetbl;
           locextfacenode_1facenode_1tbl = locextfacenode.duplicateInd({'faces', ...
                               {'faces1', 'extfaces1'}});
           locextfacenode_1facenode_1tbl = locextfacenode.duplicateInd({'nodes', ...
                               {'nodes1', 'extnodes1'}});
           % Table locextfacenode_2facenode_2tbl for Pext'
           locextfacenode_2facenode_2tbl = locextfacenodetbl;
           locextfacenode_2facenode_2tbl = locextfacenode.duplicateInd({'faces', ...
                               {'faces2', 'extfaces2'}});
           locextfacenode_2facenode_2tbl = locextfacenode.duplicateInd({'nodes', ...
                               {'nodes2', 'extnodes2'}});

           % We assemble the part corresponding to A21 = Pext*iB*div';
           
           [PextiB, locextfacenode_1facenode_2tbl] = contractTable({extsgn, ...
                               locextfacenode_1facenode_1tbl}, {iB, locface2nodetbl}, ...
                                                             {{'extfaces1', ...
                               'extnodes1'}, {'faces2', 'nodes2'}, {'faces1', ...
                               'nodes1'}});
           
           [PextiBdiv, locextfacenode_1cell_2tbl] = contractTable({PextiB, ...
                               locextfacenode_1facenode_2tbl}, {div, ...
                               loccell_2facenode_2tbl}, {{'extfaces1', ...
                               'extnodes1'}, {'cells2'}, {'faces2', ...
                               'nodes2'}});

           
           % Aggregate contribution in A21
           [~, tbl] = setupTableMapping(locextfacenode_1cell_2tbl, extfacenodetbl, ...
                                                      {{'extfaces1', 'faces'}, ...
                               {'extnodes1', 'nodes'}}); 
           locA21 = sparse(tbl.extfnind, tbl.cells2, PextiBdiv, next, nc);
           A21 = A21 - locA21;
           
           % We assemble the part corresponding to A12 = -div*iB*Pext';
           [iBPext, locfacenode_1extfacenode_2tbl] = contractTable({iB, ...
                               locface2nodetbl}, {extsgn, ...
                               locextfacenode_2facenode_2tbl}, {{'faces1', ...
                               'nodes1'}, {'extfaces2', 'extnodes2'}, {'faces2', ...
                               'nodes2'}});
           [diviBPext, loccell_1extfacenode_2tbl] = contractTable({div, ...
                               loccell_1facenode_1tbl}, {iBPext, ...
                               locfacenode_1extfacenode_2tbl}, {{'cells1'}, ...
                               {'extfaces2', 'extnodes2'}, {'faces1', ...
                               'nodes1'}});
           
           % Aggregate contribution in A12
           [~, tbl] = setupTableMapping(loccell_1extfacenode_2tbl, extfacenodetbl, ...
                                                      {{'extfaces2', 'faces'}, ...
                               {'extnodes2', 'nodes'}}); 
           locA12 = sparse(tbl.cells1, tbl.extfnind, diviBPext, nc, next);
           A12 = A12 - locA12;
           
           % Aggregate contribution in F2

           locface_1extfacenode_2tbl = projTable(locfacenode_1extfacenode_2tbl, ...
                                                 {'faces1', 'extfaces2', ...
                               'extnodes2'});
           map = setupTableMapping(locfacenode_1extfacenode_2tbl, ...
                                   locface_1extfacenode_2tbl, {'faces1', ...
                               'extfaces2', 'extnodes2'});
           iBPext = map*iBPext;            
           [~, tbl] = setupTableMapping(locface_1extfacenode_2tbl, extfacenodetbl, ...
                                                      {{'extfaces2', 'faces'}, ...
                               {'extnodes2', 'nodes'}});           
           map = setupTableMapping(locface_1extfacenode_2tbl, tbl, {'extnodes2', ...
                               'extfaces2', 'faces1'});
           iBPext = map*iBPext;
           locF2 = sparse(tbl.faces1, tbl.extfnind, iBPext, nf, next);
           
           F2 = F2 - locF2;
           
           % We assemble the part corresponding to A22 = -Pext*iB*Pext';
           [PextiB, locextfacenode_1facenode_2tbl] = contractTable({extsgn, ...
                               locextfacenode_1facenode_1tbl}, {iB, locface2nodetbl}, ...
                                                             {{'extfaces1', ...
                               'extnodes1'}, {'faces2', 'nodes2'}, {'faces1', ...
                               'nodes1'}});
           [PextiBPext, locextfacenode_1extfacenode_2tbl] = contractTable({PextiB, ...
                               locextfacenode_1facenode_2tbl}, {extsgn, ...
                               locextfacenode_2facenode_2tbl}, {{'extfaces1', ...
                               'extnodes1'}, {'extfaces2', 'extnodes2'}, ...
                               {'faces2', 'nodes2'}, });
           
           [~, tbl] = setupTableMapping(locextfacenode_1extfacenode_2tbl, ...
                                        extfacenodetbl, {{'extfaces1', 'faces'}, ...
                               {'extnodes1', 'nodes'}});
           ind1 = tbl.extfnind;
           [~, tbl] = setupTableMapping(locextfacenode_1extfacenode_2tbl, ...
                                        extfacenodetbl, {{'extfaces2', 'faces'}, ...
                               {'extnodes2', 'nodes'}});
           ind2 = tbl.extfnind;
           locA22 = sparse(ind1, ind2, PextiBPext, next, next);
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
   
   mpfastruct = struct('A'   , A, ...
                       'F'   , F, ...
                       'tbls', tbls);
end
