function [b2c_pos, b2c, locno] = invertPartition(p)
%Invert partition (cell->block mapping) to create block->cell mapping.
%
% SYNOPSIS:
%   [b2c_pos, b2c, locno] = invertPartition(p)
%
% PARAMETERS:
%   p - Partition vector.
%
% RETURNS:
%   b2c_pos - Indirection map of size [MAX(p) + 1, 1] into the 'b2c' map
%             array.  Specifically, the cells of block 'b' are stored in
%
%                  b2c(b2c_pos(b) : b2c_pos(b + 1) - 1)
%
%   b2c     - Inverse cell map.  The entries in b2c_pos(b) : b2c_pos(b+1)-1
%             correspond to the result of FIND(p == b).
%
%   locno   - A SIZE(p) array of local cell numbers.  Specifically,
%
%                 locno(c) == i
%
%             means that cell 'c' is the 'i'th cell of block p(c).
%
% NOTE:
%   This function uses SORTROWS.  Its use is therefore discouraged in
%   situations where the partition vector changes frequently, as the cost
%   of the SORTROWS call(s) may become prohibitive.
%
% EXAMPLE:
%   G = cartGrid([4, 4]);
%   p = partitionUI(G, [2, 2]);
%
%   [b2c_pos, b2c, locno] = invertPartition(p);
%
%   rot90(reshape(locno, [4, 4]))
%   plotCellData(G, locno), colorbar
%
% SEE ALSO:
%   `partitionUI`, `compressPartition`, `processPartition`, `sortrows`.

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


   c = (1 : numel(p)) .';
   a = sortrows([p, c]);
   n = accumarray(a(:,1), 1);

   b2c_pos    = cumsum([1; n]);
   b2c        = a(:,2);

   locno      = zeros(size(p));
   locno(b2c) = c - rldecode(b2c_pos(1:end-1) - 1, n);
end
