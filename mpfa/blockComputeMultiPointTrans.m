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
   if ncn > nblocks*blocksize
       blocksizes = [blocksizes; nn - nblocks*blocksize];
   end
   nblocks = numel(blocksizes);
   blockinds = cumsum([1; blocksizes]);
   
   extfaces = (G.faces.neighbors(:, 1) == 0) | (G.faces.neighbors(:, 2) == 0);
   intfaces = find(~extfaces);
   globintfacetbl.faces = intfaces;
   globintfacetbl.num   = numel(intfaces);
   
   for iblock = 1 : nblocks
       nodes = [blockinds(iblock) : (blockinds(iblock + 1) - 1)];
       [B, tbls] = blockLocalFluxMimeticAssembly(G, rock, nodes, opt);
       
       facenodetbl = tbls.facenodetbl
       [~, sz] = rlencode(facenodetbl.nodes); 
       iB   = opt.invertBlocks(B, sz);
       
       % Assemble block divergence operator, from facenode values to cell value.
       cellnodefacetbl = tbls.cellnodefacetbl;
       celltbl = tbls.celltbl;
       fno = cellnodefacetbl.faces; %alias
       cno = cellnodefacetbl.cells; %alias
       sgn = 2*(cno == G.faces.neighbors(fno, 1)) - 1;
       
       map = setupTableMapping(celltbl, cellnodefacetbl, {'cells'});
       cind = map*(celltbl.cind);
       
       div = sparse(cind, cellnodefacetbl.cnfind, sgn, celltbl.num, ...
                    cellnodefacetbl.num);
       % reduce from cell-face-node to face-node (equivalent to removing hybridization)
       op = setupTableMapping(facenodetbl, cellnodefacetbl, {'faces', 'nodes'});
       div = div*op;
   end
   
   
   %% Assemble the projection operator from facenode values to facenode values
   % on the external faces.
   fno = cellnodefacetbl.faces; %alias
   cno = cellnodefacetbl.cells; %alias
   sgn = 2*(cno == G.faces.neighbors(fno, 1)) - 1;

   extfaces = (G.faces.neighbors(:, 1) == 0) | (G.faces.neighbors(:, 2) == 0);
   faceexttbl.faces = find(extfaces);
   faceexttbl.num   = numel(faceexttbl.faces);
   [~, facenodeexttbl] = setupTableMapping(facenodetbl, faceexttbl, {'faces'});
   
   op     = setupTableMapping(cellnodefacetbl, facenodeexttbl, {'faces', 'nodes'});
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
   
   %% Assemble the system matrix operators: The degrees of freedom are the pressure
   % values at the cell center and at the external facenode.
   A11 = div*iB*div';
   A12 = -div*iB*Pext';
   A21 = Pext*iB*div';
   A22 = -Pext*iB*Pext';
   A = [[A11, A12]; [A21, A22]];

   mpfastruct = struct('iB'  , iB  , ...
                       'div' , div , ...
                       'Pext', Pext, ...
                       'F'   , F   , ...
                       'A'   , A   , ...
                       'tbls', tbls);
   
end

