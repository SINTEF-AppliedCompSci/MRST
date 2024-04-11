function stateVE = state2VE_old(state3D, Gt, fluid, res_wat, res_gas, poro3D)
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

    [S, Smax] = finescale2upscaledSat(state3D.s(:,2), res_wat, res_gas, Gt, poro3D);
    
    stateVE.s = [1-S, S];
    stateVE.sGmax = Smax;
    stateVE.pressure = finescale2upscaledPressure(state3D.pressure, Gt, fluid);
    
end
