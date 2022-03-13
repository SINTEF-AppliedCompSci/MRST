% UTILS
%
% Files
%   addCoarseCenterPoints            - Adds center points to each block. A center is typically the fine cell
%   changeCoarseCenterToWellLocation - Designates coarse nodes to cells containing wells. For coarse blocks with
%   controlVolumeRestriction         - Control volume restriction operator for a given coarsegrid partition.
%   eliminateWellEquations           - Eliminate well equations from linear system
%   formReconstructionMatrix         - Form matrix for multiscale flux reconstruction
%   getCellNeighbors                 - List of cell neighbors for a single cell.
%   getEnclosingCellsByFace          - Get cells enclosing a set of points using dot product with center-face
%   getIncomp1PhMatrix               - Get TPFA-like incompressible system matrix via incompTPFA.
%   getSystemIncompTPFA              - Extract linear system from two point solver without solving any linear
%   mapCenters                       - Map fine cells as centers of coarse cells and faces.
%   reconstructPressure              - Solve reconstruction problem for multiscale methods
%   recoverWellSolution              - Recover previously eliminated well equations from solved system
%   setCentersByWells                - Small utility to set centers of coarse blocks to one of the well cells.

%{
Copyright 2009-2022 SINTEF Digital, Mathematics & Cybernetics.

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
