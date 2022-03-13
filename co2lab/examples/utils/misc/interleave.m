function res = interleave(cellarray1, cellarray2)
% Interleave the elements of two cell arrays

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
    
  res = cell(2, max(numel(cellarray1), numel(cellarray2)));
  res(1, 1:numel(cellarray1)) = cellarray1;
  res(2, 1:numel(cellarray2)) = cellarray2;
  res = res(:);
end
