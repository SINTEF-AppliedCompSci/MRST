function binfo = identify_boundary_elements(G)
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

   sides = {'xmin', 'xmax', 'ymin', 'ymax', 'zmin', 'zmax'};
   side_faces = cell(6, 1);
   side_nodes = cell(6, 1);
   
   for i = 1:numel(sides)
      side = sides{i};
      tmp = pside([], G, side, 1);
      if isempty(tmp)
         binfo = [];
         return; % was not able to identify boundary elements
      end
      
      side_faces{i} = tmp.face;
      side_nodes{i} = ...
          unique(G.faces.nodes(mcolon(G.faces.nodePos(side_faces{i}), ...
                                      G.faces.nodePos(side_faces{i}+1)-1)));
   end
   binfo.side_faces = side_faces;
   binfo.side_nodes = side_nodes;
end
