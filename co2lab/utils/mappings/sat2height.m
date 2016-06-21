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
%   accumulateVertically, integrateVertically

%{
#COPYRIGHT#
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
