function [ncoords, ixs] = uniqueNodes(coords)
% Extract from the provided list of (non-unique) coordinates a unique set,
% and provide the corresponding indexing.
%
% SYNOPSIS:
%   function [ncoords, ixs] = uniqueNodes(coords)
%
% DESCRIPTION:
%
% PARAMETERS:
%   coords - N x D matrix of the D-dimensional coordinates for a set of
%            points.  Each row represents one point.
%
% RETURNS:
%   ncoords - M x D matrix of M points with unique nodes, resulting from
%             removing all duplicates from the 'coords' matrix.
%   ixs     - Indexing.  The point corresponding to row i in the original
%             'coords' matrix will find its coordinates at row ixs(i) in the
%             returned matrix 'ncoords'.

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

   [coords, ix] = sortrows(coords);
   
   % deciding which coordinates to keep (discarding duplicates)
   keep = [true; any(diff(coords), 2)];
   
   % Determining number of duplicates for each coord
   [~, dup] = rlencode(keep);
   
   ncoords = coords(keep,:);
   
   % Determining indexing into vector without duplicates
   ixs = zeros(size(coords, 1), 1);
   ixs(keep) = [1:sum(keep)]';
   
   diff_keep = diff(keep);
   removed_start_ixs = find(diff_keep == -1) + 1;
   removed_end_ixs   = find(diff_keep ==  1);
   if numel(removed_end_ixs) < numel(removed_start_ixs)
      removed_end_ixs = [removed_end_ixs; numel(keep)];
      assert(numel(removed_end_ixs) == numel(removed_start_ixs));
   end
   
   lengths = removed_end_ixs - removed_start_ixs + 1;
   vals = ixs(removed_start_ixs-1);
   
   ixs(~keep) = rldecode(vals, lengths);

   [~, ix] = sort(ix);
   ixs = ixs(ix);
   
end
