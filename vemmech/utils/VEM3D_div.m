function [div] = VEM3D_div(G)
% for no define all the primary operators with out units i.e without
% volum, area, length

%{ 
   Copyright 2009 - 2014 SINTEF ICT, Applied Mathematics
   %}
   assert(G.griddim == 3)
   % ncl = diff(G.cells.nodePos
   % svdiv = struct('iind', [], 'jind', [], 'matel', []); 
   % nn = sum(diff(G.faces.nodePos)) * G.griddim; 
   cells     = 1:G.cells.num; 
   ifaces    = mcolon(G.cells.facePos(cells), (G.cells.facePos(cells + 1) - 1)); 
   faces     = G.cells.faces(ifaces, 1); 
   inodes    = mcolon(G.faces.nodePos(faces), (G.faces.nodePos(faces + 1) - 1)); 
   % edges   = G.faces.edges(iedges); 
   nodes     = G.faces.nodes(inodes); 
   cellfaces = rldecode(cells', diff(G.cells.facePos)); 
   nvec      = bsxfun(@times, G.faces.normals(faces, :), ...
                      (2 * (G.faces.neighbors(faces, 1) == cellfaces) - 1)); 
   nfl       = G.faces.nodePos(faces + 1) - G.faces.nodePos(faces); 
   if(false)
      nvec = reshape(rldecode(bsxfun(@rdivide, nvec, nfl), nfl)', [], 1); % simple ???
      matel = reshape(nvec', [], 1); 
   else
      [qc, qf] = calculateQC_vec(G); 
      nvec = bsxfun(@rdivide, nvec, G.faces.areas(faces)); 
      matel = reshape(bsxfun(@times, rldecode(nvec, nfl), qf(inodes))', [], 1); 
   end
   % facenodes = rldecode(faces, diff(G.faces.nodePos)); 
   jind = mcolon(G.griddim * (nodes - 1) + 1, G.griddim * (nodes - 1) + G.griddim); 
   iind = rldecode(cellfaces, nfl * 3); 
   % assert(nn == numel(matel(:))); 
   div = sparse( iind', jind', matel(:), G.cells.num, G.nodes.num .* G.griddim);

end
