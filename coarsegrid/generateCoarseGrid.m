function cg = generateCoarseGrid(g, p, varargin)
%Form coarse grid from partition of fine-scale grid.
%
% SYNOPSIS:
%   CG = generateCoarseGrid(G, p)
%   CG = generateCoarseGrid(G, p, pf)
%
% PARAMETER:
%   G - grid_structure data structure describing fine-scale discretisation
%       of reservoir geometry.
%
%   p - Partition vector of size [G.cells.num, 1] describing the coarse
%       grid.  We assume that all coarse blocks are connected.  The
%       partition vector is often created by function partitionCartGrid
%       or function partitionUI.
%
%   pf - Partition vector on faces of 'G'. This indicator enables
%        construction of coarse grids with multiple connections between
%        coarse block pairs.  OPTIONAL.  Default value (unset) corresponds
%        to generating coarse faces defined by unique block pairs only.
%
% RETURNS:
%   CG - Coarse grid structure.  The coarse grid consists entirely of
%        topological information stored in the same topological fields as
%        in the fine-scale grid described in 'grid_structure'.
%
%        Specifically, the fields
%
%            CG.cells.num, CG.cells.facePos, CG.cells.faces
%            CG.faces.num, and CG.faces.neighbors
%
%        have the same interpretation in the coarse grid as in the
%        fine-scale grid.
%
%        The partition vector 'p' is stored within the coarse grid in a
%        separate field, 'CG.partition'.
%
% SEE ALSO:
%   coarseConnections, transposeConnections, grid_structure.

%{
Copyright 2009-2014 SINTEF ICT, Applied Mathematics.

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


   Ic = indicator(g, varargin{:});

   [conn, ind, cpos, fconn] = coarseConnections(g, p, Ic);
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

   %  3) Preserve original partition vector
   cg.partition = p;
   cg.parent    = g;
   cg.griddim   = g.griddim;
   cg.type = {'generateCoarseGrid'};
end

%--------------------------------------------------------------------------

function Ic = indicator(g, varargin)
   Ic = zeros([g.faces.num, 1]);

   if size(g.cells.faces, 2) > 1,
      % Inherit cell-face tag from FS.
      %
      ext = any(g.faces.neighbors == 0, 2);

      ix  = ext(g.cells.faces(:,1));

      Ic(g.cells.faces(ix,1)) = g.cells.faces(ix, 2);
   end

   if (nargin > 1) && ...
         (isnumeric(varargin{1}) || islogical(varargin{1})) && ...
         (size(varargin{1}, 1) == g.faces.num),

      Ic = [Ic, varargin{1}];

   end
end
