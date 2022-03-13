function cg = generateCoarseGrid(G, p, varargin)
%Form coarse grid from partition of fine-scale grid.
%
% SYNOPSIS:
%   CG = generateCoarseGrid(G, pv)
%   CG = generateCoarseGrid(G, pv, pf)
%
% PARAMETERS:
%   G -  grid_structure data structure describing fine-scale discretisation
%        of reservoir geometry.
%
%   pv - Partition vector of size [G.cells.num, 1] describing the coarse
%        grid.  We assume that all coarse blocks are connected.  The
%        partition vector is often created by function partitionCartGrid
%        or function partitionUI.
%
%   pf - Partition vector on faces of 'G'. This indicator enables
%        construction of coarse grids with multiple connections between
%        coarse block pairs.  OPTIONAL.  Default value (unset) corresponds
%        to generating coarse faces defined by unique block pairs only.
%
% RETURNS:
%   CG - Coarse grid structure. A master structure having the following
%   fields:
%    - cells --
%        A structure specifying properties for each individual block in the
%        grid.  See CELLS below for details.
%
%    - faces --
%        A structure specifying properties for each individual block
%        connections in the grid.  See FACES below for details.
%
%    - partition --
%        A copy of the partition vector 'pv'
%
%    - parent --
%        A copy of the grid structure for the parent grid
%
%    - type --
%        A cell array of strings describing the history of grid constructor
%        and modifier functions through which a particular grid structure
%        has been defined.
%
%    - griddim --
%        The dimension of the grid which in most cases will equal
%        'size(G.nodes.coords,2)'.
%
%   CELLS - Cell structure G.cells:
%    - num --
%        Number of cells in global grid.
%
%    - facePos --
%        Indirection map of size [num+1,1] into the 'cells.faces' array.
%        Specifically, the connection information of block 'i' is found in
%        the submatrix
%
%            CG.cells.faces(facePos(i) : facePos(i+1)-1, :)
%
%        The number of connections of each block may be computed using the
%        statement DIFF(facePos).
%
%    - faces --
%        A (CG.cells.facePos(end)-1)-by-2 array of global connections
%        associated with a given block.  Specifically, if
%        'cells.faces(i,1)==j', then global connection 'cells.faces(i,2)'
%        is associated with global block 'j'.
%
%        To conserve memory, only the second column is actually stored in
%        the grid structure.  The first column may be reconstructed using
%        the statement
%
%           rldecode(1 : CG.cells.num, diff(CG.cells.facePos), 2) .'
%
%        Optionally, one may append a third column to this array that
%        contains a tag that has been inherited from the parent grid.
%
%   FACES - Face structure G.faces:
%    - num --
%        Number of global connections in grid.
%
%    - connPos, fconn --
%        Packed data-array representation of coarse->fine connection
%        mapping.  Specifically, the elements
%
%            fconn(connPos(i) : connPos(i + 1) - 1)
%
%        are the connections in the parent grid (i.e., rows of the
%        neighborship definition 'G.faces.neighbors') that constitute
%        coarse-grid connection 'i'.
%
%    - neighbors --
%        A CG.faces.num-by-2 array of neighbouring information. Global
%        connection 'i' is shared by global blocks neighbors(i,1) and
%        neighbors(i,2).  One of 'neighbors(i,1)' or 'neighbors(i,2)', but
%        not both, may be zero, meaning that connection 'i' is between a
%        single block and the exterior of the grid.
%
%
%   The coarse grid consists entirely of topological information
%   stored in the same topological fields as in the fine-scale grid
%   described in 'grid_structure'. Specifically, the fields
%
%          CG.cells.num, CG.cells.facePos, CG.cells.faces
%          CG.faces.num, and CG.faces.neighbors
%
%   have the same interpretation in the coarse grid as in the
%   fine-scale grid.
%
% SEE ALSO:
%   `grid_structure`, `cellPartitionToFacePartition`, `processFacePartition`

%{
Copyright 2009-2022 SINTEF Digital, Mathematics & Cybernetics.

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


   Ic = indicator(G, varargin{:});

   [conn, ind, cpos, fconn] = coarseConnections(G, p, Ic);
   [facePos, faces]    = transposeConnections(conn);

   % Form output structure
   %  1) Cells/blocks
   cg.cells.num     = numel(facePos) - 1;
   cg.cells.facePos = facePos;
   cg.cells.faces   = [faces, ind(faces,:)];

   %  2) Interfaces
   cg.faces.num       = size(conn, 1);
   cg.faces.neighbors = conn;
   cg.faces.connPos   = cpos;
   cg.faces.fconn     = fconn;
   % Coarsen non-neighboring connections
   if isfield(G, 'nnc')
       cg.nnc = coarsenNonNeighboringConnections(G, p);
   end

   %  3) Preserve original partition vector
   cg.partition = p;
   cg.parent    = G;
   cg.griddim   = G.griddim;
   cg.type = {'generateCoarseGrid'};
end

%--------------------------------------------------------------------------

function Ic = indicator(G, varargin)
   Ic = zeros([G.faces.num, 1]);

   if size(G.cells.faces, 2) > 1,
      % Inherit cell-face tag from FS.
      %
      ext = any(G.faces.neighbors == 0, 2);

      ix  = ext(G.cells.faces(:,1));

      Ic(G.cells.faces(ix,1)) = G.cells.faces(ix, 2);
   end

   if (nargin > 1) && ...
         (isnumeric(varargin{1}) || islogical(varargin{1})) && ...
         (size(varargin{1}, 1) == G.faces.num),

      Ic = [Ic, varargin{1}];

   end
end

function nnc = coarsenNonNeighboringConnections(G, p)
    nnc = struct();
    if isfield(G.nnc, 'cells') && ~isempty(G.nnc.cells)
        % Use unique to remove merged blocks
        cells = uniqueStable(p(G.nnc.cells), 'rows');
        % Remove NNC between cells that have been placed in the same coarse
        % block.
        cells = cells(cells(:, 1) ~= cells(:, 2), :);
        nnc.cells = cells;
    end
end
