function W = reorderWellPerforationsByDepth(W, active)
%Undocumented Utility Function

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

   if nargin == 1
       active = true(numel(W), 1);
   end
   for i = 1:numel(W)
      if ~active(i)
          continue
      end
      cells   = W(i).cells;
      WI      = W(i).WI;
      dZ      = W(i).dZ;
      cstatus = W(i).cstatus;
      
      index = (1:numel(cells))';
      new_index = sortrows(horzcat(dZ, index, cells, WI, cstatus));
      W(i).dZ      = new_index(:, 1);
      W(i).cells   = new_index(:, 3);
      W(i).WI      = new_index(:, 4);
      W(i).cstatus = new_index(:, 5);
   end
end
