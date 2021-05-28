% UTILS
%
% Files
%   diagMult                         - Internal function for diagonal multiplication in AD code
%   diagProductMult                  - Undocumented Utility Function
%   double2GenericAD                 - Convert a double to GenericAD variable, using a sample GenericAD variable for dimensions
%   getSparseArguments               - Get sparse matrix indices
%   getSparseBlocks                  - Get sparse blocks
%   incrementSubset                  - Update a subset directly
%   initVariablesAD_diagonal         - Diagonal AD initializer
%   initVariablesAD_diagonalRowMajor - Diagonal AD initializer
%   initVariablesAD_oneBlock         - Initialize a set of automatic differentiation variables (single block)
%   matrixDims                       - Overloadable version of size

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
