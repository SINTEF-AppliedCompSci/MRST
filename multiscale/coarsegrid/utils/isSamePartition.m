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
Copyright 2009-2024 SINTEF Digital, Mathematics & Cybernetics.

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

      bmap = unique([reshape(p1, [], 1), reshape(p2, [], 1)], 'rows');
      b = size(bmap, 1) == size(unique(bmap(:,1)), 1);
   end
end

%--------------------------------------------------------------------------

function b = is_vector(x)
   b = (ndims(x) == 2) && (sum(size(x) == 1) >= 1);             %#ok<ISMAT>
end
