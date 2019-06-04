function mpfastruct = blockComputeNeumannMultiPointTrans(G, rock, varargin)
% Compute multi-point transmissibilities for Neumann boundary condition.
% Require incompMPFA3 solver
%
% SYNOPSIS:
%   T = computeNeumannMultiPointTrans(G, rock)
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
%   `incompMPFA3`, `mrstVerbose`.

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


   opt = struct('verbose',      mrstVerbose,   ...
                'blocksize'   , []      , ...
                'invertBlocks', 'matlab',...
                'eta', 0);

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

       
       facenodetbl = tbls.facenodetbl;
       map = setupTableMapping(intfacenodetbl, facenodetbl, {'faces', ...
                           'nodes'});
       B = map'*B*map;
       [~, sz] = rlencode(intfacenodetbl.nodes); 
       iB   = opt.invertBlocks(B, sz);
   end
   


   

   


   %% Assemble of the divergence operator, from interior facenode values to cell value.
   facenodetbl = tbls.facenodetbl;
   extfaces = (G.faces.neighbors(:, 1) == 0) | (G.faces.neighbors(:, 2) == 0);
   intfaces = find(~extfaces);
   intfacetbl.faces = intfaces;
   intfacetbl.num   = numel(intfaces);
   
   [~, intfacenodetbl] = setupTableMapping(facenodetbl, intfacetbl, {'faces'});
   % we reorder intfacenode by node so that we get the block structure for
   % the matrix B
   orderingmat = [intfacenodetbl.nodes, intfacenodetbl.faces];
   orderingmat = sortrows(orderingmat);
   intfacenodetbl.nodes = orderingmat(:, 1);
   intfacenodetbl.faces = orderingmat(:, 2);
   
   cellnodefacetbl = tbls.cellnodefacetbl;
   fno = cellnodefacetbl.faces; %alias
   cno = cellnodefacetbl.cells; %alias
   sgn = 2*(cno == G.faces.neighbors(fno, 1)) - 1;
   div = sparse(cellnodefacetbl.cells, ...
                (1 : cellnodefacetbl.num)', ...
                sgn, ...
                G.cells.num, ...
                cellnodefacetbl.num);
   % reduce from cell-face-node to face-node (equivalent to removing hybridization)
   op = setupTableMapping(intfacenodetbl, cellnodefacetbl, {'faces', 'nodes'});
   div = div*op;
   

   
   %% Invert matrix B
   % The facenode degrees of freedom, as specified by the facenodetbl table, are
   % ordered by nodes first (see implementation below). It means in particular
   % that the matrix B is, by construction, block diagonal.

   if opt.verbose
       fprintf('Computing inverse mixed innerproduct...\n');
       t0 = tic();   
   end

   map = setupTableMapping(intfacenodetbl, facenodetbl, {'faces', ...
                       'nodes'});
   
   B = map'*B*map;
   [~, sz] = rlencode(intfacenodetbl.nodes); 
   iB   = opt.invertBlocks(B, sz);
   
   if opt.verbose
       t0 = toc(t0);   
       fprintf('... computing inverse mixed innerproduct done in %g sec\n', t0);
   end
   
   
   %% Assemble the flux operator: From pressure values at the cell center and
   % at the external facenode, compute the fluxes at the faces
   F = iB*div';
   map = setupTableMapping(intfacenodetbl, intfacetbl, {'faces'});
   F = map*F;
   
   %% Assemble the system matrix operaror: The degrees of freedom are the pressure
   % values at the cell center and at the external facenode.
   A = div*iB*div';

   mpfastruct = struct('iB'  , iB  , ...
                       'div' , div , ...
                       'F'   , F   , ...
                       'A'   , A   , ...
                       'tbls', tbls);
   
end

