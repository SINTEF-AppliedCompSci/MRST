% Files
%   addFaces.m                   - Add faces F from grid structure
%   computeTransportSourceTerm.m - Compute source terms for transport
%   copyFaces.m                  - Copy face data for faces F from grid structure
%   implicitTransport.m          - Implicit single point upwind transport solver for two-phase flow.
%   insertInPackedData.m         - Insert c into row r of packed array (data, pos)
%   matrixBlocksFromSparse.m     - Extract block-diagonal matrix elements from sparse matrix
%   newtonRaphson2ph.m           - Solve non-linear equation F(s)=0 using Newton-Raphson method.
%   removeFaces.m                - Remove faces F from grid structure
%   removeFromPackedData.m       - Remove all copies of specific elements from packed data structure

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
