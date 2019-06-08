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
%   `incompMPFA2`, `mrstVerbose`.

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

   opt = struct('verbose'     , mrstVerbose, ...
                'blocksize'   , []      , ...
                'invertBlocks', 'matlab', ...
                'eta'         , 0);

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
   globintfacetbl.faces = intfaces;
   globintfacetbl.num   = numel(intfaces);

   cellfacetbl.cells = rldecode((1 : nc)', diff(G.cells.facePos));
   cellfacetbl.faces = G.cells.faces(:, 1);
   cellfacetbl.num   = numel(cellfacetbl.cells);

   facenodetbl.faces = rldecode((1 : nf)', diff(G.faces.nodePos));
   facenodetbl.nodes = G.faces.nodes;
   facenodetbl.num = numel(facenodetbl.faces);

   [~, cellfacenodetbl] = setupTableMapping(cellfacetbl, facenodetbl, ...
                                                         {'faces'});
   A = sparse(nc, nc);

   for iblock = 1 : nblocks

       nodes = [blockinds(iblock) : (blockinds(iblock + 1) - 1)]';
       [B, tbls] = blockLocalFluxMimeticAssembly(G, rock, nodes, opt);

       locfacenodetbl     = tbls.facenodetbl;
       locface2nodetbl    = tbls.face2nodetbl;
       loccellfacenodetbl = tbls.cellfacenodetbl;
       loccellnodetbl     = tbls.cellnodetbl;

       % Assembly of B
       Bmat = sparse(locface2nodetbl.fnind1, locface2nodetbl.fnind2, B, ...
                     locfacenodetbl.num, locfacenodetbl.num);
       [~, sz] = rlencode(locfacenodetbl.nodes);
       iBmat   = opt.invertBlocks(Bmat, sz);
       clear locmattbl
       locmattbl.fnind = (1 : locfacenodetbl.num)';
       locmattbl.num = locfacenodetbl.num;
       [~, locmattbl] = setupTableMapping(locmattbl, locmattbl, [], 'duplicate', ...
                                                     {{'fnind', {'fnind1', ...
                           'fnind2'}}});
       locmattbl = sortTable(locmattbl, {'fnind2', 'fnind1'});
       iB = iBmat(:);
       map = setupTableMapping(locmattbl, locface2nodetbl, {'fnind1', 'fnind2'});
       iB = map*iB;
       clear map;

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

       % Table loccell_1face_1nodetbl for div mapping from facenode to cell
       loccell_1face_1nodetbl = replacefield(loccellfacenodetbl, 'faces', ...
                                           'faces1');
       loccell_1face_1nodetbl = replacefield(loccell_1face_1nodetbl, 'cells', ...
                                             'cells1');
       loccell_1face_1nodefds = {'cells1', 'faces1', 'nodes'};

       % Table loccell_2face_2nodetbl for div' mapping from cell to facenode
       loccell_2face_2nodetbl = replacefield(loccellfacenodetbl, 'faces', ...
                                           'faces2');
       loccell_2face_2nodetbl = replacefield(loccell_2face_2nodetbl, 'cells', ...
                                             'cells2');
       loccell_2face_2nodefds = {'cells2', 'faces2', 'nodes'};

       [~, prodmattbl] = setupTableMapping(loccell_1face_1nodetbl, locface2nodetbl, ...
                                                       {'faces1', 'nodes'});
       map1 = setupTableMapping(loccell_1face_1nodetbl, prodmattbl, {'faces1', ...
                           'cells1', 'nodes'});
       map2 = setupTableMapping(locface2nodetbl, prodmattbl, {'faces1', ...
                           'nodes'});
       locA = (map1*div).*(map2*iB);

       % Table loccell_1face_2nodetbl for div*iB mapping (-> locA) from facenode to cell
       loccell_1face_2nodetbl = projTable(prodmattbl, {'cells1', 'faces2', ...
                           'nodes'});
       loccell_1face_2nodefds = {'cells1', 'faces2', 'nodes'};

       reducemap = setupTableMapping(prodmattbl, loccell_1face_2nodetbl, {'cells1', ...
                           'faces2', 'nodes'});
       locA = reducemap*locA; % mapping in loccell_1face_1nodetbl

       % figure
       % tbl1 = projTable(loccell_1face_2nodetbl, {'cells1'});
       % tbl1 = addLocInd(tbl1, 'cind');
       % tbl2 = projTable(loccell_1face_2nodetbl, {'faces2', 'nodes'});
       % tbl2 = addLocInd(tbl2, 'fnind');
       % map1 = setupTableMapping(tbl1, loccell_1face_2nodetbl, {'cells1'});
       % ind1 = map1*tbl1.cind;
       % map2 = setupTableMapping(tbl2, loccell_1face_2nodetbl, {'faces2', 'nodes'});
       % ind2 = map2*tbl2.fnind;
       % locAmat = sparse(ind1, ind2, locA, tbl1.num, tbl2.num);
       % spy(locAmat);

       [~, prodmattbl] = setupTableMapping(loccell_2face_2nodetbl, ...
                                        loccell_1face_2nodetbl, {'faces2', ...
                           'nodes'});
       map1 = setupTableMapping(loccell_1face_2nodetbl, prodmattbl, ...
                                              {'cells1', 'faces2', 'nodes'});
       map2 = setupTableMapping(loccell_2face_2nodetbl, prodmattbl, ...
                                              {'cells2', 'faces2', 'nodes'});
       locA = (map1*locA).*(map2*div);
       clear map1 map2

       % Table loccell2tbl for div*iB*div' (-> locA) from cell to cell
       loccell2tbl = projTable(prodmattbl, {'cells1', 'cells2'});
       loccell2fds = {'cells1', 'cells2'};
       reducemap = setupTableMapping(prodmattbl, loccell2tbl, {'cells1', ...
                           'cells2'});
       locA = reducemap*locA;
       clear reducemap

       % Assemble sparse matrix
       locA = sparse(loccell2tbl.cells1, loccell2tbl.cells2, locA, nc, nc);

       A = A + locA;

   end

   tbls = struct('facenodetbl', facenodetbl, 'cellfacenodetbl', ...
                 cellfacenodetbl);

   [bcstructs, tbls] = boundaryStructures(G, rock, tbls, 'eta', opt.eta);

   F    = bcstructs.F;
   Aold = bcstructs.A;
   A12  = bcstructs.A12;
   A21  = bcstructs.A21;
   A22  = bcstructs.A22;

   A11  = A;
   A = [[A11, A12]; [A21, A22]];
   
   mpfastruct = struct('F'   , F  , ...
                       'A'   , A  , ...
                       'A11' , A11, ...
                       'tbls', tbls);
end

function [bcstructs, tbls] = boundaryStructures(G, rock, tbls, varargin)

    opt = struct('verbose'     , mrstVerbose, ...
                 'invertBlocks', 'matlab', ...
                 'eta'         , 0);

    opt = merge_options(opt, varargin{:});
    opt.invertBlocks = blockInverter(opt);

    [B, rtbls] = robustComputeLocalFluxMimetic(G, rock, opt);

    %% Invert matrix B
    % The facenode degrees of freedom, as specified by the facenodetbl table, are
    % ordered by nodes first (see implementation below). It means in particular
    % that the matrix B is, by construction, block diagonal.
    rfacenodetbl = rtbls.facenodetbl;
    [~, sz] = rlencode(rfacenodetbl.nodes);
    iB   = opt.invertBlocks(B, sz);

    facenodetbl = tbls.facenodetbl;
    map = setupTableMapping(facenodetbl, rfacenodetbl, {'faces', 'nodes'});
    iB = map'*iB*map;

    %% Assemble of the divergence operator, from facenode values to cell value.
    cellfacenodetbl = tbls.cellfacenodetbl;
    fno = cellfacenodetbl.faces; %alias
    cno = cellfacenodetbl.cells; %alias
    sgn = 2*(cno == G.faces.neighbors(fno, 1)) - 1;
    div = sparse(cellfacenodetbl.cells, ...
                 (1 : cellfacenodetbl.num)', ...
                 sgn, ...
                 G.cells.num, ...
                 cellfacenodetbl.num);
    % reduce from cell-face-node to face-node (equivalent to removing hybridization)
    op = setupTableMapping(facenodetbl, cellfacenodetbl, {'faces', 'nodes'});
    div = div*op;

    %% Assemble the projection operator from facenode values to facenode values
    % on the external faces.
    fno = cellfacenodetbl.faces; %alias
    cno = cellfacenodetbl.cells; %alias
    sgn = 2*(cno == G.faces.neighbors(fno, 1)) - 1;

    extfaces = (G.faces.neighbors(:, 1) == 0) | (G.faces.neighbors(:, 2) == 0);
    faceexttbl.faces = find(extfaces);
    faceexttbl.num   = numel(faceexttbl.faces);
    [~, facenodeexttbl] = setupTableMapping(facenodetbl, faceexttbl, {'faces'});

    op     = setupTableMapping(cellfacenodetbl, facenodeexttbl, {'faces', 'nodes'});
    fn_sgn = op*sgn;
    map = setupTableMapping(facenodetbl, facenodeexttbl, {'faces', 'nodes'});
    nfne = facenodeexttbl.num;
    Pext = sparse(1 : nfne, 1 : nfne, fn_sgn, nfne, nfne)*map;

    tbls.facenodeexttbl = facenodeexttbl;

    %% Assemble the flux operator: From pressure values at the cell center and
    % at the external facenode, compute the fluxes at the faces
    F1 = iB*div';
    F2 = - iB*Pext';
    F  = [F1, F2];
    facetbl.faces = (1 : G.faces.num)';
    facetbl.num   = G.faces.num;
    map = setupTableMapping(facenodetbl, facetbl, {'faces'});
    F = map*F;

    %% Assemble the system matrix operaror: The degrees of freedom are the pressure
    % values at the cell center and at the external facenode.
    A11 = div*iB*div';
    A12 = -div*iB*Pext';
    A21 = Pext*iB*div';
    A22 = -Pext*iB*Pext';
    A = [[A11, A12]; [A21, A22]];

    bcstructs = struct('iB'  , iB  , ...
                       'div' , div , ...
                       'Pext', Pext, ...
                       'F'   , F   , ...
                       'A'   , A   , ...
                       'A11' , A11 , ...
                       'A12' , A12 , ...
                       'A21' , A21 , ...
                       'A22' , A22);
end
