function [qc, qf, qcvol] = calculateQC_vec(G)
% calculate the q used informulation so it is in G.cells.nodes format

%{ 
   Copyright 2009 - 2014 SINTEF ICT, Applied Mathematics
%} 
   assert(G.griddim == 3); 
   % qf = nan(size(G.faces.nodes)); 
   % qc = zeros(numel(G.cells.nodes), 3); 
   % triangles and regular tetrahedrons is q = 1 / n; 
   % make face node - cellnode mapping.
   
   cellno       = rldecode([1:G.cells.num]', diff(G.cells.nodePos)); % #ok
   A            = sparse(cellno, G.cells.nodes, 1:numel(G.cells.nodes)); 

   faces        = [1:G.faces.num]'; 

   inodes       = mcolon(G.faces.nodePos(faces), G.faces.nodePos(faces + 1) - 1); 
   nodes        = G.faces.nodes(inodes); 
   nlc          = diff(G.faces.nodePos); 
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
   % EL         = sqrt(sum(E.^2, 2)); 
   NE           = bsxfun(@cross, E, rldecode(N, nlc)); 
   % NE         = bsxfun(@rdivide, NE, EL); 
   qf           = zeros(numel(inodes), 1); 
   dxbb         = rldecode(xb - xbb, nlc); 
   qf(ind1_e)   = sum((NE(ind2_e, :) + NE(ind1_e, :)) .* dxbb(ind1_e, :), 2); 
   qf(ind1)     = sum((NE(ind2, :) + NE(ind1, :)) .* dxbb(ind1, :), 2); 
   qf           = rldecode(G.faces.areas(faces) ./ nlc, nlc) + 0.5 * qf; 

   % qf(inodes) = ql; 
   %% 
   num_cn = numel(G.cells.nodes); 
   qcvol = zeros(num_cn,1);
   qc = zeros(numel(G.cells.nodes), G.griddim); 
   normals = bsxfun(@rdivide, G.faces.normals, G.faces.areas); 
   for j = 1:2
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
      rcord = fcord-ccord;%vector from cells centroid to centroid/point on face

      % try to map to local numbering of nodes in cell
      lind  = sub2ind(size(A), cells, nodes(ind)); 
      cn    = reshape(full(A(lind)), [], 1); 
      %
      pyrvol = sum(rcord.*qcl,2);
      qcvol = qcvol+accumarray(cn, pyrvol, [num_cn, 1])/3;
      for i = 1:G.griddim
         qc(:, i) = qc(:, i) + accumarray(cn, qcl(:, i), [num_cn, 1]); 
      end
      
   end


   %{
   % do cell assemble iin loop at the moment
   for face = 1:G.faces.num
      inodes = G.faces.nodePos(face):G.faces.nodePos(face + 1) - 1; 
      nodes  = G.faces.nodes(inodes); 
      nn     = numel(nodes); 
      ql     = qf(inodes); 
      % 
      N = G.faces.normals(face, :) ./ G.faces.areas(face); 
      ncell = G.faces.neighbors(face, :); 
      if(ncell(1) == 0)
         ncell(1) = ncell(2); 
         N = -N; 
         ncell(2) = 0; 
      end
      cn1 = nan(nn, 1); 
      cn2 = nan(nn, 1); 
      for i = 1:nn
         cn1(i) = A(ncell(1), nodes(i)); 
         if(ncell(2)>0)
            cn2(i) = A(ncell(2), nodes(i)); 
         end
      end
      qN = bsxfun(@times, N, ql); 
      qc(cn1(:), :) = qc(cn1(:), :) + qN; 
    if(ncell(2)>0)    
        qc(cn2(:), :) = qc(cn2(:), :)-qN;
    end
end
end
%}
