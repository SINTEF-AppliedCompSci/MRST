function [cno, nno, hfno, fno, subfno, subcno] = createSubcellTopology(g)
%
% Mostly cut and paste from subfunction createMappings in
% computeMultiPointTrans; however not the final output
%
   % Create mapping from sub-half-face to cell, node, face, half-face and
   % sub-face
   
   % Cell numbers corresponding to g.cells.faces
   cellno   = rldecode(1:g.cells.num, diff(g.cells.facePos), 2) .';
   % In which column in g.faces.neighbors is the cell
   col      = 1+(cellno == g.faces.neighbors(g.cells.faces(:,1), 2));
   nhfaces  = g.cells.facePos(end)-1;
   % (g.faces.num) x 2 matrix, telling where in g.cells.faces each face
   % occurs (as here and there in the columns)
   hfaces   = accumarray([g.cells.faces(:,1), col], 1:nhfaces);
   % Duplicate the information, so that there is one entry for each node
   % the face is connected to
   hfaces   = rldecode(hfaces, diff(g.faces.nodePos));

   % Duplicate the neighbors of each face, one for each node the face is
   % connected to
   cells    =  rldecode(g.faces.neighbors, diff(g.faces.nodePos));
   % Nodes of the faces, one time for each side
   nodes    =  repmat(g.faces.nodes, [2,1]);
   % Face index, corresponding to the nodes
   faces    =  repmat(rldecode(1:g.faces.num, diff(g.faces.nodePos),2)', [2,1]);
   % Define an index for subfaces (note that this is the only thing new
   % here, the other quantities are simply reconstructed from the topology)
   subfaces =  repmat((1:size(g.faces.nodes,1))', [2,1]);
   % Zero cell means we are outside the boundary
   i        =  cells~=0;
   w        =  [cells(i), nodes(i), hfaces(i), faces(i), subfaces(i)];
   w        =  double(sortrows(w));


   cno     = w(:,1);
   nno     = w(:,2);
   hfno    = w(:,3);
   fno     = w(:,4);
   subfno  = w(:,5);
   
   % Associate an index with each subcell (i.e. combination of cell and
   % node)
   [~,~,subcno] = unique([cno,nno],'rows');
   subhfno = (1:numel(cno))';
end