function sG = getGasSaturationFromHeight(T, t, B, b, h)
    % Internal function: Reconstruct gas saturation from VE assumption
    % T is top of column for each cell, t is top of cells to be
    % reconstructed, with b/B having the same interpretation for the
    % respective cell bottoms. h is the height of the gas plume in each
    % cell.

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

    d_local = t + h;
    sG = (d_local - T)./(B-T);
    sG = max(sG, 0); % Gas < 0 means that the gas-liquid interface is actually above the cell top
    sG = min(sG, 1); % Gas > 1 means that the gas-liquid interface is below the cell bottom
end
