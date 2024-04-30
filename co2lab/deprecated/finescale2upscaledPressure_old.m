function P = finescale2upscaledPressure_old(p, Gt, fluid)
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
 
    % index of uppermost cell in each column
    cells_upper = Gt.columns.cells(Gt.cells.columnPos(1:end-1));
    
    p_upper = p(cells_upper);
    
    dz_upper = Gt.columns.dz(cells_upper);
    
    rho = fluid.rhoWS .* fluid.bW(p_upper);
    
    P = p_upper - rho .* dz_upper/2 * norm(gravity());
    
end
