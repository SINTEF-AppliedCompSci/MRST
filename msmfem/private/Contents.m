% Files
%   assignBasisFuncs.m    - Update reservoir coarse system basis functions following 'evalBasisFunc'.
%   checkBoundaryBasis.m  - Check structural integrity of coarse system versus boundary conditions.
%   evalBasisFunc.m       - Compute multiscale basis functions for selected faces.
%   evalBasisFuncGlobal.m - Compute multiscale basis functions for selected faces.
%   evalWellBasis.m       - Compute multiscale basis functions for single well.
%   splitBlocks.m         - Split a set of coarse blocks into sub-blocks.
%   splitBlocksForWells.m - Split coarse blocks near wells.

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
