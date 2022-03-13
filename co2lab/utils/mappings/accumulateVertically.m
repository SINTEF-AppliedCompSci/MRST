function q = accumulateVertically(q, p)
%Accumulate quantity vertically per column.
%
% SYNOPSIS:
%   cq = accumulateVertically(q, pos)
%
% PARAMETERS:
%   q   - Sampled quantity values.  One scalar value for each cell in each
%         column.
%
%   pos - Column indirection array.  Specifically,
%
%             pos(i) : pos(i + 1) - 1
%
%         are the indices in 'q' which correspond to column 'i'.
%
% RETURNS:
%   cq - Accumulated 'q' values.  Corresponds to CUMSUM(q), but reset per
%        column such that ALL(cq(pos(1 : end-1)) == 0).
%
% SEE ALSO:
%   `topSurfaceGrid`, `cumsum`.

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

% $Date: 2012-01-30 11:39:51 +0100 (Mon, 30 Jan 2012) $
% $Revision: 9019 $

   ix = 1 : numel(p) - 1;
   q  = arrayfun(@(i) reshape(cumsum(q(p(i) : p(i+1) - 1)), [], 1), ...
                 ix, 'UniformOutput', false);
   q  = vertcat(q{:});
end
