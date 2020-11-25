function G = createBisectedTetrahedralGrid(Nx,alternate)
%Undocumented Utility Function

%{
Copyright 2020 University of Bergen and SINTEF Digital, Mathematics & Cybernetics.

This file is part of the MPSA-W module for the MATLAB Reservoir Simulation Toolbox (MRST).

The MPSA-W module is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

The MPSA-W module is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with the MPSA-W module.  If not, see <http://www.gnu.org/licenses/>.
%}


    [X,Y,Z] = meshgrid(linspace(0,1,Nx(1)+1), linspace(0,1,Nx(2)+1),linspace(0,1,Nx(3)+1));
    p = [X(:), Y(:), Z(:)];

    G = mcomputeGeometry(tetrahedralGrid(p));
end
