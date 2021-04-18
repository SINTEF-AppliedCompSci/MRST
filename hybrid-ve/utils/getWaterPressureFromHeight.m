function pW_local = getWaterPressureFromHeight(pt, t, B, b, h, pW, g, rhow, rhog)
    % Internal function: Reconstruct water pressure from VE assumption
    % T is top of column for each cell, t is top of cells to be
    % reconstructed, with b/B having the same interpretation for the
    % respective cell bottoms. h is the height of the gas plume in each
    % cell.
    % pW is the water pressure at the bottom of the cell, g the gravity
    % magnitude and rhow/rhog the water and gas densities.

%{
Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.

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

    dz_below = (b - pt);

    h_below = max(t + h - pt, 0);
    pW_local = pW - ((dz_below - h_below).*rhow + h_below.*rhog).*g;
end
