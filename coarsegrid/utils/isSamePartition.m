function b = isSamePartition(p1, p2)
%Check if two partition vectors represent the same grid partition
%
% SYNOPSIS:
%   b = isSamePartition(p1, p2)
%
% PARAMETERS:
%   p1, p2 - Partition vectors.
%
% RETURNS:
%   b - Whether or not partitions p1 and p2 are exactly equal or, failing
%       that, whether or not one is a renumbering/permutation of the other.
%       Logical scalar.
%
% NOTE:
%   This function uses sortrows in all but the most trivial cases.
%
% EXAMPLE:
%   p1 = partitionCartGrid([9, 9, 9], [3, 3, 3]);
%   i  = randperm(max(p1));
%   p2 = 10 * i(p1);              % Note: Holes in partition vector.
%   b  = isSamePartition(p1, p2)

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

   % Trivial case: p1 == p2
   b = isequal(p1, p2);

   if ~b && isnumeric(p1) && isnumeric(p2) && ...
         is_vector(p1) && is_vector(p2) && ...
         (numel(p1) == numel(p2))
      % p1 ~= p2, but both numeric vectors and both have the same number of
      % elements.  Check if one is uncompressed permutation of the other.
      %
      % Note: Extracting block-to-cell mapping requires SORTROWS.

      [b2c_p1, b2c_1] = block_cells(p1);
      [b2c_p2, b2c_2] = block_cells(p2);

      if sum(diff(b2c_p1) > 0) == sum(diff(b2c_p2) > 0)
         % Same number of active blocks.  Possibly same partition.
         %
         % Check if fine-scale cells partition equally in both vectors.
         %
         % Note: By construction, cell indices appear in sorted order in
         % the b2c array.
         cells = @(b, pos, b2c) ...
            reshape(b2c(pos(b) : (pos(b + 1) - 1)), [], 1);

         c1 = @(b1) cells(b1, b2c_p1, b2c_1);
         c2 = @(b2) cells(b2, b2c_p2, b2c_2);
         eq = @(c)  isequal(c1(p1(c)), c2(p2(c)));

         % Loop representative (first) cells from each block of p1 and see
         % if the blocks containing those cells have the same contents in
         % both p1 and p2.
         b = all(arrayfun(eq, b2c_1(b2c_p1(diff(b2c_p1) > 0))));
      else
         % Different number of active blocks so p1 and p2 do definitely not
         % represent the same cell partitions.  Note: This assignment is,
         % strictly speaking, a no-op because 'b' is already FALSE from the
         % earlier ISEQUAL check.
         b = false;
      end
   end
end

%--------------------------------------------------------------------------

function b = is_vector(x)
   b = (ndims(x) == 2) && (sum(size(x) == 1) >= 1);             %#ok<ISMAT>
end

%--------------------------------------------------------------------------

function [pos, b2c] = block_cells(p)
   [pos, b2c] = invertPartition(reshape(p, [], 1));
end
