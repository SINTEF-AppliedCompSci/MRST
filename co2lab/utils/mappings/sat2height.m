function h = sat2height(s, g, rock)
%Convert from saturation to height
%
% SYNOPSIS:
%   h = sat2height(s, G, rock)
%
% PARAMETERS:
%   s  - Saturation - one value for each cell in the underlying 3D model.
%        Corresponds to state.s for the 3D problem.
%
%   G  - A top-surface grid as defined by function 'topSurfaceGrid'.
%
%   rock - 3D rock data containing valid filed 'poro'.
%
% RETURNS:
%   h - Plume thickness.  One scalar value for each column in the top-surface 
%       grid.
%
% SEE ALSO:
%   `accumulateVertically`, `integrateVertically`

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

% pore height in each cell in 3D model
phCell = rock.poro(g.columns.cells).*s(g.columns.cells).*g.columns.dz;

% compute cumulative quantities
ph = cumulativePoreHeight(g,rock);

colNo = rldecode((1:g.cells.num)', diff(g.cells.columnPos));

% pore height in each column in 2D model
phCol = accumarray(colNo, phCell);

h = invertVerticalFunction(phCol, g, g.columns.z,  ph);
end
