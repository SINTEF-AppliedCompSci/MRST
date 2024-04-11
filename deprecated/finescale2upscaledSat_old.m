function [S, Smax] = finescale2upscaledSat_old(sg, res_wat, res_gas, Gt, poro)
%Undocumented Utility Function

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

    cix = Gt.columns.cells;
    pvol = poro(cix) .* Gt.parent.cells.volumes(cix);
    gvol = pvol .* sg(cix);
    
    cmap = rldecode((1:Gt.cells.num)', diff(Gt.cells.columnPos));
    col_pvol = accumarray(cmap, pvol);
    col_gvol = accumarray(cmap, gvol);
    
    S = col_gvol ./ col_pvol; % volume of gas in a column as a fraction of
                              % the total porevolume of the column
    
    % to determine Smax, we fill all cells with gas saturation >= res_gas 
    % completely with CO2, and compute S again
    threshold = 0.9; % allow residual saturation to be slightly lower than res_gas
    sgmax = sg;
    % full saturation in all cells with saturation about residual threshold
    sgmax(sgmax > threshold * res_gas) = 1 - res_wat; 
    
    gmaxvol = pvol .* sgmax(cix);
    col_gmaxvol = accumarray(cmap, gmaxvol);
    Smax = col_gmaxvol ./ col_pvol;
    
end

