function [div] = VEM2D_div(G, varargin)
% for no define all the primary operators with out units i.e without
% volum, area, length


%{ 
   Copyright 2009 - 2014 SINTEF ICT, Applied Mathematics
   %} 
   opt = struct('addface_dof', false); 
   opt = merge_options(opt, varargin{:}); 
   div = VEM2d_loc_int(G, opt); 
end

% ---------------------------------------------------------------------------- 

function div = VEM2d_loc_int(G, opt)
   qf     = calculateQF_vec(G); 
   qf     = reshape(qf', [], 1); 
   dofs   = mcolon(G.griddim * (G.cells.nodes - 1) + 1, ...
                   G.griddim * (G.cells.nodes - 1) + G.griddim); 
   ndofs  = G.nodes.num * G.griddim; 
   cellno = rldecode([1:G.cells.num]', diff(G.cells.nodePos) * G.griddim); 
   
   if(opt.addface_dof)
      
      ncf        = diff(G.cells.facePos); 
      fcellno    = rldecode([1:G.cells.num]', ncf); 
      faces      = G.cells.faces(:, 1); 
      % nfacedof = numel(faces); 
      fasgn      = (2 * (G.faces.neighbors(faces, 2) == fcellno) - 1) .* ...
                         G.faces.areas(faces); 
      cellno     = [cellno; fcellno]; dofs = [dofs'; ndofs + faces]; 
      qf         = [qf; fasgn]; 
      ndofs      = ndofs + G.faces.num; 
   
   end    
   
   div = sparse(cellno, dofs, qf, G.cells.num, ndofs); 
   
end

% ---------------------------------------------------------------------------- 

function div = VEM2d_loc(G)
% calculate the q used informulation so it is in G.cells.nodes format
   qf = zeros(size(G.cells.faces, 1), 2); 
   % qc = zeros(numel(G.cells.nodes), 3); 
   % triangles and regular tetrahedrons is q = 1 / n; 
   % make face node - cellnode mapping.
   cellno = rldecode([1:G.cells.num]', diff(G.cells.nodePos)); % #ok

   % A = sparse(cellno, G.cells.nodes, 1:numel(G.cells.nodes)); 
   for numblocks = 1:1% G.cells.num
      cells = 1:G.cells.num; 
      inodes1 = mcolon(G.cells.nodePos(cells), G.cells.nodePos(cells + 1) - 2); 
      inodes2 = mcolon(G.cells.nodePos(cells) + 1, G.cells.nodePos(cells + 1) - 1); 
      ifaces = mcolon(G.cells.facePos(cells), G.cells.facePos(cells + 1) - 1); 
      faces = G.cells.faces(ifaces, 1); 
      sign = 2 * (G.faces.neighbors(faces, 1) == cellno) - 1; 
      N = bsxfun(@times, G.faces.normals(faces', :), sign); 
      qf(inodes1, :) = qf(inodes1, :) + N(inodes1, :); 
      qf(inodes2, :) = qf(inodes2, :) + N(inodes1, :); 
      qf(G.cells.nodePos(cells), :) = ... 
          qf(G.cells.nodePos(cells), :) + ...
          N(G.cells.nodePos(cells + 1) - 1, :); 
      qf(G.cells.nodePos(cells + 1) - 1, :) = ...
          qf(G.cells.nodePos(cells + 1) - 1, :) + ...
          N(G.cells.nodePos(cells + 1) - 1, :); 
      qf = qf / 2; 
      inodes = mcolon(G.cells.nodePos(cells), G.cells.nodePos(cells + 1) - 1); 
      nodes = G.cells.nodes(inodes); 
      dofs = mcolon(2 * nodes - 1, 2 * nodes); 
      cellnum = rldecode(cells', diff(G.cells.nodePos) * 2); 
      qf = qf'; 
      div = sparse(cellnum, dofs', qf(:), G.cells.num, 2 * G.nodes.num); 
   end
end
