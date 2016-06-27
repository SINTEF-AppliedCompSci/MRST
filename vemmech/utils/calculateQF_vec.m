function [qf, qf_vol] = calculateQF_vec(G)
% calculate the q used informulation so it is in G.cells.nodes format

%{ 
   Copyright 2009 - 2014 SINTEF ICT, Applied Mathematics
   %}
   assert(G.griddim == 2); 
   
   % qf has one entry per node per cell
   qf = zeros(size(G.cells.nodes, 1), 2); 
   %qf = zeros(size(G.cells.faces, 1), 2); @@ previously used

   % qc = zeros(numel(G.cells.nodes), 3); 
   % triangles and regular tetrahedrons is q = 1 / n; 
   % make face node - cellnode mapping.

   % Cell indices repeated according to how many nodes (or faces) it has.  An
   % implicit assumption is that the number of nodes and the number of faces
   % for a cell is always equal (should hold in 2D).
   cellno = rldecode([1:G.cells.num]', diff(G.cells.nodePos)); % #ok

   % A = sparse(cellno, G.cells.nodes, 1:numel(G.cells.nodes)); 
   for numblocks = 1:1 % G.cells.num
      cells   = 1:G.cells.num;
      lcells = rldecode(cells',diff(G.cells.nodePos)');
      
      % For each cell, indices of the first to the second-to-last node 
      % (indirection level 2)
      inodes1 = mcolon(G.cells.nodePos(cells), G.cells.nodePos(cells + 1) - 2)'; 
      
      % For each cell, indices of the second to the last node (indirection level 2)
      inodes2 = mcolon(G.cells.nodePos(cells) + 1, G.cells.nodePos(cells + 1) - 1)'; 

      % For each cell, indices of each face (indirection level 2)
      ifaces  = mcolon(G.cells.facePos(cells), G.cells.facePos(cells + 1) - 1)'; 
      
      % For each cell, indices of each face (indirection level 1)
      faces   = G.cells.faces(ifaces, 1); 
      
      % Orienting normals ('N') so that they always point out of the current
      % cell and into the neighbor cell
      sign    = 2 * (G.faces.neighbors(faces, 1) == cellno) - 1; 
      N       = bsxfun(@times, G.faces.normals(faces', :), sign); 

      % For each cell node, add up the (scaled) normals of the two adjacent faces and
      % divide by two.
      relvec = G.faces.centroids(faces,:)-G.cells.centroids(lcells,:);
      tetvols = sum(N.*relvec,2);
      qf_vol = zeros(numel(G.cells.nodes),1);
      qf_vol(inodes1) = qf_vol(inodes1, :) + tetvols(inodes1); 
      qf_vol(inodes2) = qf_vol(inodes2) +  tetvols(inodes1);
      qf_vol(G.cells.nodePos(cells)) = qf_vol(G.cells.nodePos(cells)) + ...
                                              tetvols(G.cells.nodePos(cells+1)-1);
      qf_vol(G.cells.nodePos(cells + 1) - 1) = qf_vol(G.cells.nodePos(cells + 1) -1) + ...
                                              tetvols(G.cells.nodePos(cells+1)-1);                                    
      qf_vol = qf_vol/4;
      
      qf(inodes1, :) = qf(inodes1, :) + N(inodes1, :); 
      qf(inodes2, :) = qf(inodes2, :) + N(inodes1, :);
      qf(G.cells.nodePos(cells), :) = qf(G.cells.nodePos(cells), :) + ...
                                      N(G.cells.nodePos(cells + 1) - 1, :); 
      qf(G.cells.nodePos(cells + 1) - 1, :) = qf(G.cells.nodePos(cells + 1) -1,:) + ...
                                              N(G.cells.nodePos(cells+1)-1,:);
      qf = qf / 2;
      % The two components for each node (line) in qf thus represent the
      % integral of the basis function in each coordinate over the non-zero
      % faces for that basis function.
end

