function ph = cumulativePoreHeight(g, rock)
%Compute cumulative pore height for each column.
%
% SYNOPSIS:
%   ph = cumulativePoreHeight(G, rock)
%
% DESCRIPTION:
%   Define the pore height of a cell c as dz(c)*phi(c).
%
% PARAMETERS:
%   G    - A top-surface grid as defined by function 'topSurfaceGrid'.
%
%   rock - Rock data structure.  Must contain valid field 'rock.poro'.
%
% RETURNS:
%   pI - A G.cells.num-by-1 array of doubles such that
%
%           pv(G.cells.columnPos(i) : G.cells.columnPos(i+1)-1))
%
%        contains the cumulative pore height, measured from top-level z==0,
%        at the bottom interface of each 3D cell in column 'i'.
%
% NOTE:
%   This function is to the pore height what function cumulativeHeight is
%   to the incremental depth, G.columns.dz.
%
% SEE ALSO:
%   `topSurfaceGrid`, `cumulativeHeight`, `accumulateVertically`.

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

   assert (isfield(g        , 'columns'  ));
   assert (isfield(g.columns, 'dz'       ));
   assert (isfield(g.cells  , 'columnPos'));
   assert (isfield(rock     , 'poro'     ));

   poro = rock.poro(g.columns.cells);
   ph   = accumulateVertically(poro .* g.columns.dz, g.cells.columnPos);
end
