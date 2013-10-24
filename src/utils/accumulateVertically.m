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
%   topSurfaceGrid, cumsum.

%{
#COPYRIGHT#
%}

% $Date: 2012-01-30 11:39:51 +0100 (Mon, 30 Jan 2012) $
% $Revision: 9019 $

   ix = 1 : numel(p) - 1;
   q  = arrayfun(@(i) reshape(cumsum(q(p(i) : p(i+1) - 1)), [], 1), ...
                 ix, 'UniformOutput', false);
   q  = vertcat(q{:});
end
