function z = cumulativeHeight(g)
%Compute cumulative height for each column.
%
% SYNOPSIS:
%   z = cumulativeHeight(G)
%
% PARAMETERS:
%   G - A top-surface grid as defined by function 'topSurfaceGrid'.
%
% RETURNS:
%   z - A G.cells.num-by-1 array of doubles such that
%
%           z(G.columns.pos(i) : G.columns.pos(i+1)-1))
%
%       contains the cumulative height, measured from top-level z==0, of
%       the bottom interface of each 3D cell in column 'i'.  Specifically,
%       when
%
%           p = G.columns.pos(1 : end-1),
%
%       then
%
%           z(p) == G.columns.dz(p)
%
%       modulo arithmetic error (round-off).
%
% SEE ALSO:
%   topSurfaceGrid, accumulateVertically.

%{
#COPYRIGHT#
%}

% $Date: 2012-01-30 11:39:51 +0100 (Mon, 30 Jan 2012) $
% $Revision: 9019 $

   assert (isfield(g, 'columns'));
   assert (isfield(g.columns, 'dz'));
   assert (isfield(g.cells, 'columnPos'));

   z = accumulateVertically(g.columns.dz, g.cells.columnPos);
end
