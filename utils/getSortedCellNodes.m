function cn = getSortedCellNodes(G)
% Construct n x 2 table of cell edges with edges oriented the same
% direction around the cell boundary.

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

   cellNo    = rldecode(1:G.cells.num, diff(G.cells.facePos), 2).';
   edges     = reshape(G.faces.nodes, 2, [])';
   cellEdges = edges(G.cells.faces(:,1),:);
   ind       = G.faces.neighbors(G.cells.faces(:,1), 1) ~= cellNo;
   cellEdges(ind, :) = cellEdges(ind, [2,1]);

   % Check for existing sortedCellNodes and do basic consistency check
   if isfield(G.cells, 'sortedCellNodes')
       cn = G.cells.sortedCellNodes;
       assert(size(cn, 1) == size(cellEdges, 1), ...
           ['sortedCellNodes is malformatted. Remove field' ...
           'G.cells.sortedCellNodes or call getSortedCellNodes again.'])
       return
   end

   % Sort edges in each cell:
   for c = 1 : G.cells.num
      ind = G.cells.facePos(c) : G.cells.facePos(c + 1) - 1;
      cellEdges(ind, :) = sortEdges(cellEdges(ind,:));
   end
   cn = reshape(cellEdges(:,1), 1, [])';
end

%--------------------------------------------------------------------------

function edges = sortEdges(edges)
   % Assume edges vectors are oriented in the same direction around cell.
   % Sort edges such that they are back-to-back.
   % Then cellNodes are edges(:,1).

   for i = 1 : size(edges, 1) - 1
      for j = i + 1 : size(edges,1)

         % Check if j follows edges(i)
         if any(edges(i, 2) == edges(j, :), 2)
            if edges(i,2) == edges(j,2), edges(j,:) = edges(j,[2,1]);end
            % Add j to end of list by swapping edges i+1 and j
            tmp = edges(i+1,:);
            edges(i+1, :) = edges(j,:);
            edges(j,   :) = tmp;
            break
         end

         % Check if j precedes edges(1)
         if any(edges(1, 1) == edges(j, :), 2)
            if edges(1,1) == edges(j,1), edges(j,:) = edges(j,[2,1]);end
            % Add j to front of list
            edges = edges([j,1:j-1,j+1:end],:);
            break
         end

      end
   end
end
