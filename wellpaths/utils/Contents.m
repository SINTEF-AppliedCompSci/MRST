% UTILS
%
% Files
%   combineWellPaths   - Combine multiple simple paths into a full tree
%   findWellPathCells  - Convert well path to the intersected cells
%   getWellFromPath    - Convert well path to MRST well.
%   makeSingleWellpath - Create well path from points (and optional connectivity and active flag)
%   plotWellPath       - Plot a well path
%   refineSpline       - Refine a curve to higher resolution using spline interpolation

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
