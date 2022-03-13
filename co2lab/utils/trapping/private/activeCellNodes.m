function cnodes = activeCellNodes(Gt)
% Return a (4 x m)-sized matrix where m is the number of active cells in Gt.  Each
% column holds the indices of the 4 nodes that are corners of that cell.  
% The algorithm presupposes that each cell has exactly 4 corner nodes.
    
    cnode_full_mat = cellNodes(Gt);
    
    % only keep the information we want
    cnodes = reshape(cnode_full_mat(:, 3), 4, []);
       
end

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
