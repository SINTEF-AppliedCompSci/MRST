function binfo = identify_boundary_elements(G)

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
