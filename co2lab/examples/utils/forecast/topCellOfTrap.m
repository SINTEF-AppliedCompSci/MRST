function tops = topCellOfTrap(Gt, ta)
% explicitly get top cells of each trap region
% (in case ta.tops is not indexed properly according to trap number)

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

    num_traps = numel(ta.trap_z);
    tops = [];
    for i = 1:num_traps
        tcells = find(ta.traps == i);
        %tz = zeros(Gt.cells.num,1);
        tz = NaN(Gt.cells.num,1);
        tz(tcells) = Gt.cells.z(tcells);
        top = find(tz == min(tz)); % top at the min elevation of trap i
        tops(i,1) = top(1);  % the top elevation may occur at more than one cell, so we only
                             % need to take one cell to get the top
        clear tcells
    end
end
