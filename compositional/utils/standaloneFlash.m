function [L, x, y, Z_L, Z_V] = standaloneFlash(p, T, z, EOSModel)
% Utility for flashing without explicitly forming a state
    solver = NonLinearSolver();
    solver.maxTimestepCuts = 0;
    solver.maxIterations = 200;
    solver.continueOnFailure = true;
    solver.errorOnFailure = false;
    
    state = struct();

    state.pressure = p;
    state.T = T;
    state.components = z;
    
    state = solver.solveTimestep(state, 1, EOSModel);
    
    L = state.L;
    x = state.x;
    y = state.y;
    Z_L = state.Z_L;
    Z_V = state.Z_V;
end

%{
Copyright 2009-2016 SINTEF ICT, Applied Mathematics.

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
