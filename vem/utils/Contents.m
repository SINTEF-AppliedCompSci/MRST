% UTILS
%
% Files
%   addBCVEM            - Add boundary condition to (new or existing) BC using function handels,
%   computeVEMGeometry  - Computes VEM geometry of MRST grid G.
%   conserveFlux        - Postprocess nonconservative flux field.
%   polyDim             - Computes the dimension of the space of polynomials of degree k or less in
%   polygonInt          - Integrates the function f over each cell in cells of grid G, using a
%   polygonInt3D        - Integrates the function f over each face in faces of 3D grid G, using a
%   polyhedronInt       - Integrates the function f over each cell in cells of grid G, using a
%   squeezeBlockDiag    - Squeezes a block diagonal matrix in which each block has the same number
%   tetrahedronQuadRule - Returns quadrature rule for the reference terahedron with vertices V.
%   triangleQuadRule    - Returns quadrature rule for the reference triangle with vertices (0,0),

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
