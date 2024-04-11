function p3D = pVE_to_p3D(Gt, pVE, fluid)
% Convert from upscaled pressure to fine-scale pressure.
%   
% SYNOPSIS:
%   p3D = pVE_to_p3D(Gt, pVE, fluid)
% 
% PARAMETERS:
%   Gt    - the top surface grid on which the upscaled pressure field is defined 
%   pVE   - the upscaled (2D) pressure field, representing the pressure at the 
%           caprock level
%   fluid - the VE fluid object used during simulation (needed here to
%           compute brine density).
%
% RETURNS:
%   p3D - reconstructed 3D pressure field, defined on the original 3D grid
%         domain (found in Gt.parent).

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

    gravity_state = norm(gravity) > 0; % check if gravity is on
    gravity on; % turn on gravity temporarily (in case it was turned off)
    
    rhoW = fluid.rhoW(pVE); % water density at top surface
    
    % number of cells in each pillar
    cells_per_column = diff(Gt.cells.columnPos);

    % for each cell in the 3D grid, identify the corresponding cell in the 
    % top surface grid:

    % list of topcells, repeated for each cell in a column
    topcells = rldecode((1:Gt.cells.num)', cells_per_column);
    
    % attribute correct topcell to each 3D grid cell
    topcells = topcells(Gt.columns.cells); 
    
    top_z_val = Gt.cells.z(topcells); % z position of pVE for the cells
                                      % in each column
    
    G = Gt.parent;
    
    z_offsets  = G.cells.centroids(:,3) - top_z_val; % z offset from VE pressure
                                                     % reference to individual cell depths
    
    p3D = pVE(topcells) + rhoW(topcells) .* norm(gravity()) .* z_offsets;
    
    gravity(gravity_state); % set gravity state back to what it was
end
