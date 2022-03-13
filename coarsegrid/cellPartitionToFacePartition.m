function pf = cellPartitionToFacePartition(g, p, varargin)
%Construct partition of all grid faces from cell partition.
%
% SYNOPSIS:
%   pf = cellPartitionToFacePartition(g, p)
%   pf = cellPartitionToFacePartition(g, p, 'pn1', pv1, ...)
%
% DESCRIPTION:
%   Define unique partition ID for each pair of cell partition IDs
%   occurring in a partition vector, and construct a partitioning of all
%   faces in a grid.
%
% PARAMETERS:
%   g       - Grid structure as described by grid_structure.
%
%   p       - Cell partition
%
%   'pn'/pv - List of 'key'/value pairs for supplying optional parameters.
%             The supported options are
%                AllBoundaryFaces -
%                  Flag indicating whether or not every fine-scale boundary
%                  face should be assigned a separate, unique partition ID
%                  number to force individual inclusion in the coarse grid.
%
%                  Logical.  Default value: AllBoundaryFaces=false (do not
%                  assign separate partition IDs to each fine-grid boundary
%                  face).
%
% RETURNS:
%   pf      - Face partition.  One non-negative integer per face in 'g'.
%             Fine-scale faces that are not on the boundary between coarse
%             blocks are assigned a partition ID number of zero.
%
% SEE ALSO:
%   `processFacePartition`

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


   opt = struct('AllBoundaryFaces', false);
   opt = merge_options(opt, varargin{:});

   % Form face partition pf from cell partition p.
   p       = [0;p];
   [tmp,k] = sortrows(sort(p(g.faces.neighbors+1), 2));
   [n, n]  = rlencode(tmp);                                            %#ok
   pf      = zeros(g.faces.num, 1);
   pf(k)   = rldecode(1:numel(n), n, 2)';

   % Add all fine-faces on domain boundary
   if opt.AllBoundaryFaces,
      ix     = any(g.faces.neighbors==0, 2);
      pf(ix) = max(pf)+(1:sum(ix))';
   end
end
