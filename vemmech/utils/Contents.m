% UTILS
%
% Files
%   addFluidContribMechVEM - SYNOPSIS:
%   C2D                    - Compute the matrix of normalized strain energies (D) from the elasticity
%   calculateQC            - SYNOPSIS:
%   calculateQF            - SYNOPSIS:
%   complex3DGrid          - SYNOPSIS:
%   Enu2C                  - SYNOPSIS:
%   ENu2LMu_3D             - SYNOPSIS:
%   exploreSquareGrid      - Explore the different types of grids that can be set up using the function squareGrid
%   lincompTPFA            - Solve weakly compressible flow problem (fluxes/pressures) using TPFA method.
%   linearDisplacement     - SYNOPSIS:
%   squareGrid             - SYNOPSIS:
%   squareLayersTest       - SYNOPSIS:
%   squareTest             - Different test cases for linear elasticity on square domains
%   VEM2D_div              - Discrete div operator for the virtual element method in 2D
%   VEM3D_div              - Discrete div operator for the virtual element method in 3D
%   VEM_div                - Discrete div operator for the virtual element method

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
