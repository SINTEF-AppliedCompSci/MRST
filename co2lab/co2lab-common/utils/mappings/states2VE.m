function statesVE = states2VE(states3D, Gt, fluid, poro3D)
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

    sGmax = zeros(Gt.cells.num, 1);
    statesVE = cell(numel(states3D), 1);
    
    for i = 1:numel(states3D)
    
        state = states3D{i};
        sg = finescale2upscaledSat(state.s(:,2), Gt, poro3D);
        p  = finescale2upscaledPressure(state.pressure, Gt, fluid);
        
        sGmax = max(sGmax, sg);
        
        sVE.pressure = p;
        sVE.s = [1-sg, sg];
        sVE.sGmax = sGmax;
        statesVE{i} = sVE;
    end
    
end

% ----------------------------------------------------------------------------
function S = finescale2upscaledSat(sg, Gt, poro)

    cix = Gt.columns.cells;
    pvol = poro(cix) .* Gt.parent.cells.volumes(cix);
    gvol = pvol .* sg(cix);
    
    cmap = rldecode((1:Gt.cells.num)', diff(Gt.cells.columnPos));
    col_pvol = accumarray(cmap, pvol);
    col_gvol = accumarray(cmap, gvol);
    
    S = col_gvol ./ col_pvol; % volume of gas in a column as a fraction of
                              % the total porevolume of the column
end

% ----------------------------------------------------------------------------
function P = finescale2upscaledPressure(p, Gt, fluid)
 
    % index of uppermost cell in each column
    cells_upper = Gt.columns.cells(Gt.cells.columnPos(1:end-1));
    
    p_upper = p(cells_upper);
    
    dz_upper = Gt.columns.dz(cells_upper);
    
    rho = fluid.rhoWS .* fluid.bW(p_upper);
    
    P = p_upper - rho .* dz_upper/2 * norm(gravity());
    
end
