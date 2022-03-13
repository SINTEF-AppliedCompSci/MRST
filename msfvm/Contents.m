% MSFVM
%
% Files
%   createPermutationMatrix - Create permutation matrices for the coarse orderings in the MsFV
%   find_clusters           - Undocumented Utility Function
%   findCenter              - indices of blocks in the current cell
%   mldivide_update         - A,b linear system
%   partitionUIdual         - Create coarse dual partition/grid
%   plotDual                - Plot an implicitly defined dual grid for the multiscale finite volume method
%   reportError             - Simple error helper for MsFV examples
%   solveMSFV_TPFA_Incomp   - Solve incompressible flow problem (flux/pressures) using a multiscale

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
