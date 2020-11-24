function G = gridForConvTest(Nx,gridType)
%Undocumented Utility Function

%{
Copyright 2009-2020 SINTEF Digital, Mathematics & Cybernetics.

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

    %{
      Grid type:
      1: Cartesian
      2: Triangles by alternating bisection of triangles
      3: Equilateral triangles
      4: Triangles by uniform bisection
      5: Tetrehedral grid  
    %}

    Nd = numel(Nx);

    switch gridType
      case 1
        G = computeGeometry(cartGrid(Nx,ones(1,Nd)));
      case 2
        G = createBisectedTriangleGrid(Nx,1);
      case 3
        G = createEquilateralTriangleGrid(Nx);
      case 4
        G = createBisectedTriangleGrid(Nx,0);
      case 5
        assert(Nd == 3)
        G = createBisectedTetrahedralGrid(Nx);
      otherwise
        error('gridType not recognized');
    end
end
