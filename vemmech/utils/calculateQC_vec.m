function [qc, qf, qcvol] = calculateQC_vec(G)
%
% SYNOPSIS:
%   function [qc, qf, qcvol] = calculateQC_vec(G)
%
% DESCRIPTION:  Calculate elementary integrals that are used to assemble the
% stiffness matrix for the 3D case. The precise definitions can be found in [Gain et al].
%
% PARAMETERS:
%   G - Grid structure
%
% RETURNS:
%   qc    - Elementary assembly integrals : One (3D) vector value in each
%                                           cell, see (74) in [Gain et al].  
%   qf    - Elementary assembly integrals : One scalar value for each face in each cells,
%                                           corresponds to (98) in [Gain et al].
%   qcvol - Elementary assembly integrals : One scalar value for each node in each cells,
%                                           gives weights to compute the L^2
%                                           projections, see VEM_linElast.m
%
% EXAMPLE:
%
% SEE ALSO:
%



% calculate the q used informulation so it is in G.cells.nodes format

%{ 
   Copyright 2009 - 2014 SINTEF ICT, Applied Mathematics
%} 
   assert(G.griddim == 3); 
   
   cellno       = rldecode([1:G.cells.num]', diff(G.cells.nodePos)); % #ok
   A            = sparse(cellno, G.cells.nodes, 1:numel(G.cells.nodes)); 

   faces        = [1 : G.faces.num]'; 

   inodes       = mcolon(G.faces.nodePos(faces), G.faces.nodePos(faces + 1) - 1); 
   nodes        = G.faces.nodes(inodes); 
   nlc          = diff(G.faces.nodePos);  % Number of nodes per cell
   facenode     = rldecode(faces, nlc); 
   xn           = G.nodes.coords(nodes, :); 
   xbb          = zeros(G.faces.num, G.griddim); 

   for i  = 1:G.griddim
      xbb(:, i) = accumarray(facenode, G.nodes.coords(nodes, i)); 
   end

   xbb          = bsxfun(@rdivide, xbb, nlc); 
   xb           = G.faces.centroids(faces, :); 
   N            = bsxfun(@rdivide, G.faces.normals(faces, :), G.faces.areas(faces)); 
   E            = zeros(size(xn)); 
   ind1         = mcolon(G.faces.nodePos(faces) + 1, G.faces.nodePos(faces + 1) - 1); 
   ind2         = mcolon(G.faces.nodePos(faces), G.faces.nodePos(faces + 1) - 2); 
   E(ind2, :)   = -xn(ind2, :) + xn(ind1, :); 
   ind2_e       = G.faces.nodePos(faces + 1) - 1; 
   ind1_e       = G.faces.nodePos(faces); 
   E(ind2_e, :) = -xn(ind2_e, :) + xn(ind1_e, :); 
   NE           = bsxfun(@cross, E, rldecode(N, nlc)); 
   qf           = zeros(numel(inodes), 1); 
   dxbb         = rldecode(xb - xbb, nlc); 
   qf(ind1_e)   = sum((NE(ind2_e, :) + NE(ind1_e, :)) .* dxbb(ind1_e, :), 2); 
   qf(ind1)     = sum((NE(ind2, :) + NE(ind1, :)) .* dxbb(ind1, :), 2); 
   qf           = rldecode(G.faces.areas(faces) ./ nlc, nlc) + 0.5 * qf; 

   num_cn = numel(G.cells.nodes); 
   qcvol = zeros(num_cn,1);
   qc = zeros(numel(G.cells.nodes), G.griddim); 
   normals = bsxfun(@rdivide, G.faces.normals, G.faces.areas); 
   for j = 1 : 2
      cells = rldecode(G.faces.neighbors(:, j), nlc); 
      qcl   = bsxfun(@times, rldecode(normals, nlc), qf);
      fcord = rldecode(G.faces.centroids, nlc);
      if(j == 2)
         qcl = -qcl; 
      end
      ind   = cells>0; 
      cells = cells(ind); 
      qcl   = qcl(ind, :);
      fcord = fcord(ind,:);
      ccord = G.cells.centroids(cells,:);
      rcord = fcord-ccord; % vector from cells centroid to centroid/point on face

      % map to local numbering of nodes in cell
      lind  = sub2ind(size(A), cells, nodes(ind)); 
      cn    = reshape(full(A(lind)), [], 1); 
      %
      pyrvol = sum(rcord.*qcl, 2);
      qcvol = qcvol + accumarray(cn, pyrvol, [num_cn, 1])/3;
      for i = 1 : G.griddim
         qc(:, i) = qc(:, i) + accumarray(cn, qcl(:, i), [num_cn, 1]); 
      end
   end
end
