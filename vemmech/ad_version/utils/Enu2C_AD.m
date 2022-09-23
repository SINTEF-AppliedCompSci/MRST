function C = Enu2C_AD(E, nu, G, ixnames)
%Undocumented Utility Function

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

   assert(G.griddim == 3) % @@ only 3D supported for now
   if ~exist('ixnames', 'var')
      ixnames = {'i', 'j', 'c'}; % default index names
   end
   
   
   % scalar factor
   T1 = SparseMultiArray(E ./ (1 + nu) ./ (1 - 2 * nu));
   
   % upper-left block
   T2 = SparseMultiArray(nu * ones(9, 1), ...
                     vertcat(repmat([1,2,3], 1, 3), ...
                             [ones(1,3), 2*ones(1,3), 3*ones(1,3)])', ...
                     ixnames(1:2));
   
   T3 = SparseMultiArray((1 - 2 * nu) * ones(3, 1), ...
                     [1, 1; 2, 2; 3, 3], ixnames(1:2));
   
   % lower-right block
   T4 = SparseMultiArray((1 - 2 * nu)/2 * ones(3, 1), ...
                     [4, 4; 5, 5; 6, 6], ixnames(1:2));
   
   % tensor representing cell-space
   Tcells = SparseMultiArray([], (1:G.cells.num)', ixnames(3));
   
   % resulting C tensor
   
   C = (T1 * (T2 + T3 + T4)) * Tcells;
end
