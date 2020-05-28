function [S, operators] = VEM_assemble_AD(G, E, nu, varargin)

   opt = merge_options(struct('extra', []), varargin{:});

   %% Compute non-AD part
   % compute matrices that do not depend on the (potentially AD) elastic
   % moduli, but only on geometric or toplogical information. 
   % These matrices include: 'I', 'volmap', 'Wc', 'Nc', 'PP' and 'assemb'
   % If 'extra' is given as an optional parameter, these matrices will not be
   % recomputed, only forwarded.
   non_AD_part = compute_geometry_dependent_matrices(G, opt.extra);
   
   %% Compute AD part
   % Compute tensors that depend on the grid and elastic moduli, but not on
   % any of the matrices computed above.
   AD_part = compute_moduli_dependent_tensors(G, E, nu);
   
   %% Assemble the final matrices
   % Compute the tensors/matrices that depend on both the non-AD and the AD
   % components computed above
   [S, operators] = final_assembly(non_AD_part, AD_part);
   
end

% ----------------------------------------------------------------------------
function res = compute_geometry_dependent_matrices(G, pcomp)
   
   % If precomputed values are provided, return these and skip the rest of the 
   if ~isempty(pcomp)
      [res.I, res.volmap, res.WC, res.NC, res.PP, res.assemb] = ...
          deal(pcomp.I, pcomp.volmap, pcomp.WC, pcomp.NC, pcomp.PP, pcomp.assemb);
      return;
   end
   
   % If we got here, precomputed values were not provided. We must compute
   % them from the grid.

   % first, compute some basic indexing stuff
   cells = 1:G.cells.num;
   inodes = mcolon(G.cells.nodePos(cells), (G.cells.nodePos(cells + 1)-1))';
   nodes = G.cells.nodes(inodes);
   linodes = (1:numel(inodes))';
   nldofs = G.griddim * numel(inodes); % number of linear degs. of freedom
   nlc = G.cells.nodePos(cells + 1) - G.cells.nodePos(cells); % # of nodes per cell
   lcellnum = rldecode((1:numel(cells))', nlc);

   if G.griddim == 2
      nlin = 3; % dimension ov Voigt matrix (and of local linear space)
   else
      nlin = 6;
   end
   
   % compute some basic geometric quantities
   BB = NaN(numel(cells), 1);  % cell centers
   for i = 1:G.griddim
      BB(:, i) = accumarray(lcellnum, G.nodes.coords(nodes, i), [numel(cells), 1]); 
   end
   XB = G.nodes.coords(nodes, :) - BB(lcellnum, :); % coords relative to center
   
   % define identity matrix of the full system size
   I = sparse(1:nldofs, 1:nldofs, 1, nldofs, nldofs);
   
   % compute cell volumes associated with cell-node pairs ('volmap')
   vol = G.cells.volumes(cells); % cell volumes
   volfac = rldecode(vol, nlc * G.griddim);
   volmap = sparse(1:nldofs, 1:nldofs, volfac, nldofs, nldofs); 
   
   % construct NC and NR, which give the nodal representations of the local linear
   % transformations
   nn = numel(nodes);
   [zz, z] = deal(zeros(nn, 2), zeros(nn, 1));
   [NC, NR] = deal(zeros(G.griddim * nn, nlin));
      
   if(G.griddim == 3)
      NC(:, 1:3) = [reshape([XB(:, 1), zz]', [], 1), ...
                    reshape([z, XB(:, 2), z]', [], 1), ...
                    reshape([zz, XB(:, 3)]', [], 1)];
      perm = {[2, 1, 3], [1, 3, 2], [3, 2, 1]};
      nrv = {[1, -1, 0], [0, 1, -1], [-1, 0, 1]};
      nrc = {[1, 1, 0], [0, 1, 1], [1, 0, 1]};
   else
      NC(:, 1:2) = [reshape([XB(:, 1), z]', [], 1), ...
                    reshape([z, XB(:, 2)]', [], 1)];
      perm = {[2, 1]};
      nrv = {[1, -1]};
      nrc = {[1, 1]};
   end
   NR(:, 1:G.griddim) = repmat(eye(G.griddim), nn, 1);
   for i = 1:numel(perm)
      tmp = XB(:, perm{i});
      NR(:, G.griddim + i) = reshape(bsxfun(@times, tmp, nrv{i})', [], 1);
      NC(:, G.griddim + i) = reshape(bsxfun(@times, tmp, nrc{i})', [], 1);
   end
        
   % construct WC and WR, which compute the local linear tranformations from given
   % nodal displacements
   [WC, WR] = deal(zeros(G.griddim * nn, nlin));
   
   nlcl = rldecode(nlc, nlc);
   cellnum = rldecode(cells, nlc);
   
   if (G.griddim == 3)
      qc_all = calculateQC(G);
      % XX is 2*q_i in [Gain et al : doi:10.1016/j.cma.2014.05.005] eqs 77
      XX = bsxfun(@rdivide, qc_all(inodes, :), G.cells.volumes(cellnum));
      WC(:, 1:3) = [reshape([XX(:, 1), zz]', [], 1), ...
                    reshape([z, XX(:, 2), z]', [], 1), ...
                    reshape([zz, XX(:, 3)]', [], 1)];
      XX = XX / 2; % now XX is q_i
      WR(:, 1:3) = [reshape([1./nlcl, zz]', [], 1), ...
                    reshape([z, 1./nlcl, z]', [], 1), ...
                    reshape([zz, 1./nlcl]', [], 1)];
   else
      qf_all = calculateQF(G);
      XX = bsxfun(@rdivide, qf_all(inodes, :), G.cells.volumes(cellnum));
      WC(:, 1:2) = [reshape([XX(:, 1), z]', [], 1), ...
                    reshape([z, XX(:, 2)]', [], 1)];
      XX = XX / 2;
      WR(:, 1:2) = [reshape([1./nlcl, z]', [], 1), ...
                    reshape([z, 1./nlcl]', [], 1)];
   end
      
   for i = 1:numel(perm)
      tmp = XX(:, perm{i});
      WR(:, G.griddim + i) = reshape(bsxfun(@times, tmp, nrv{i})', [], 1);
      WC(:, G.griddim + i) = reshape(bsxfun(@times, tmp, nrc{i})', [], 1);
   end
   
   % convert the matrices NC, NR, WC andWR to full-sized block matrices
   [ib, jb] = blockDiagIndex(nlin * ones(size(nlc)), G.griddim * nlc);
   mat_sz = numel(nodes) * G.griddim;

   WR = sparse(ib, jb, reshape(WR', [], 1), nlin*numel(cells), mat_sz)';
   WC = sparse(ib, jb, reshape(WC', [], 1), nlin*numel(cells), mat_sz)';
   NR = sparse(ib, jb, reshape(NR', [], 1), nlin*numel(cells), mat_sz)';
   NC = sparse(ib, jb, reshape(NC', [], 1), nlin*numel(cells), mat_sz)';
   
   % define the complete projection matrices 
   PR = NR * WR'; % rigid body rotation
   PC = NC * WC'; % other linear displacements
   PP = PR + PC; % projection to linear displacement
   
   % assembly matrix ('assemb')
   dofsg = mcolon(G.griddim * (nodes - 1) + 1, G.griddim * (nodes - 1) + G.griddim);
   dofsl = mcolon(G.griddim * (linodes - 1) + 1, G.griddim * (linodes - 1) + G.griddim);
   assemb = sparse(dofsg, dofsl, 1, G.griddim * G.nodes.num, numel(linodes) * G.griddim);
   
   % organize the return values in a return structure
   [res.I, res.volmap, res.WC, res.NC, res.PP, res.assemb] = ...
       deal(I, volmap, WC, NC, PP, assemb);
end

% ----------------------------------------------------------------------------

function res = compute_moduli_dependent_tensors(G, E, nu)
   
   D = ;
   SE = ;
   
   res.D = D;
   res.SE = SE;
end

% ----------------------------------------------------------------------------

function [S, operators] = final_assembly(nad, ad)
  
  
end
