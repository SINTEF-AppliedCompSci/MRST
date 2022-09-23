function [S, operators] = VEM_assemble_AD(G, E, nu,  varargin)
% Assemble the virtual element stiffness matrix and intermediate operators.
% The function is similar to `VEM_assemble`, but material properties 'E' and
% 'nu' are allowed to be ADI-variables, and the matrix and operators are
% returned in the form of SparseMultiArrays.
%
% SYNOPSIS:
%   function [S, operators] = VEM_assemble_AD(G, C, varargin)
%
% DESCRIPTION:
%   Compute the full VEM stiffness matrix; optionally returns intermediate
%   operators
% 
%   Assembly based on Paulino's paper in 3D. The notations follow Paulino's
%   paper.
%
%   Vectorized version. 
% 
%   S = |E| W_c D W_c^T + (I - P_P)^T s (I - P_P)
%
% PARAMETERS:
%   G        - Grid (2D or 3D) 
%   C        - matrix where each row represents the elasticity tensor for 
%              the grid cell corresponding to that row
%   varargin - options are:
%              'blocksize' - size of blocks (# of cells) for vectorized
%              calc.
%              'alpha_scaling' - scaling of the stabilization term.  Default
%              is 1.
%              'extra' - a previously computed version of the 'operators'
%              structure can be passed back in to this function for re-use,
%              to avoid having to re-assemble the parts that do not depend on
%              material properties E and nu.
%
% RETURNS:
%   S         - full VEM system matrix, including rows/columns for
%   Dirichlet nodes. operators - will contain the fields:
%               - D      - matrix of normalized strain energies for each
%                          cell
%               - WC     - matrix for projection of basis functions onto
%                          space of constant strain modes 'c' 
%               - assemb -
%
% EXAMPLE:
%
% SEE ALSO: `VEM_linElast`

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
   
   opt = merge_options(struct('extra'        , [], ...
                              'alpha_scaling', 1 ), ...
                       varargin{:});

   %% Compute non-AD part
   % compute matrices that do not depend on the (potentially AD) elastic
   % moduli, but only on geometric or toplogical information. 
   % These matrices include: 'I', 'volmap', 'Wc', 'Nc', 'PP' and 'assemb'
   % If 'extra' is given as an optional parameter, these matrices will not be
   % recomputed, only forwarded.    
   % If 'salt' is not empty, these matrices will include a modified WC
   % system with only compressive modes
   fprintf('-- non_AD_part.\n');
   non_AD_part = compute_geometry_dependent_matrices(G, opt.extra);
   
   %% Compute AD part
   % Compute the tensor that depend on the grid and elastic moduli, but not on
   % any of the matrices computed above.
   fprintf('-- AD part.\n');
   D = compute_moduli_dependent_tensor(G, E, nu);
   
   %% Assemble the final matrices
   % Compute the tensors/matrices that depend on both the non-AD and the AD
   % components computed above
   fprintf('-- final assembly\n');
   [S, operators] = final_assembly(G, non_AD_part, D, opt.alpha_scaling);
   fprintf('-- assembly finished.\n');
end

% ----------------------------------------------------------------------------
function res = compute_geometry_dependent_matrices(G,precomp)
   
   % If precomputed values are provided, return these and skip the rest of the 
   if ~isempty(precomp)
      [res.I, res.volmap, res.WC, res.NC, res.PP, res.assemb] = ...
          deal(precomp.I, precomp.volmap, precomp.WC, ...
               precomp.NC, precomp.PP, precomp.assemb);
      return;
   end
   
   % If we got here, precomputed values were not provided. We must compute
   % them from the grid.

   % first, compute some basic indexing stuff
   cells = (1:G.cells.num)';
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
   fac = accumarray(lcellnum, 1, [numel(cells), 1]);
   BB = bsxfun(@rdivide, BB, fac);
   
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

function D = compute_moduli_dependent_tensor(G, E, nu)
   
   if G.griddim==2
      nufac = [-1,  1, 0; 
                1, -1, 0; 
                0,  0, -1];% nu-factor
      cdiag = diag([1,1, 1/2]); % constant diagonal
      fac = [1,1,2]';
   else % G.griddim==3
      nufac = [-1,  1,  1,  0,  0,  0;
                1, -1,  1,  0,  0,  0;
                1,  1, -1,  0,  0,  0;
                0,  0,  0, -1,  0,  0;
                0,  0,  0,  0, -1,  0;
                0,  0,  0,  0,  0, -1];
      cdiag = diag([1, 1, 1, 1/2, 1/2, 1/2]);
      fac = [1,1,1,2,2,2]';
   end

   % C is the elasticity tensor
   C = (SparseMultiArray(nufac, {'i', 'j'}) * SparseMultiArray(nu, {'cell'}) + ...
        ( SparseMultiArray(cdiag, {'i', 'j'}) * SparseMultiArray([], (1:G.cells.num)', {'cell'}))) ^ ...
       SparseMultiArray(E ./ (1 + nu) ./ (1 - 2 * nu), {'cell'});
   
   % D is the material property tensor used in the VEM formulation
   D = C ^ SparseMultiArray(fac, 'i') ^ SparseMultiArray(fac, 'j');

end

% ----------------------------------------------------------------------------

function [S, op] = final_assembly(G, op, D, alpha_scaling)

   % initial indexing stuff
   if G.griddim == 2
      nlin = 3; % dimension ov Voigt matrix (and of local linear space)
   else
      nlin = 6;
   end
   cells = (1:G.cells.num)';
   nlc = G.cells.nodePos(cells + 1) - G.cells.nodePos(cells); % # of nodes per cell
   inodes = mcolon(G.cells.nodePos(cells), (G.cells.nodePos(cells + 1)-1))';
   nldofs = G.griddim * numel(inodes); % number of linear degs. of freedom   
   
   % compute trace of D
   fprintf('----Compute trace of D\n');
   trD = D.contract('i', 'j').asVector({'cell'});
   
   % compute trace of Nc'Nc
   fprintf('----Compute trace of Nc''Nc');
   DNC = diag(op.NC' * op.NC);
   trDNC = sum(reshape(DNC, nlin, []), 1)';  
   
   % cannot use rldecode since alpha could be AD.  Use sparse matrix instead
   %   alpha = rldecode((trD .* G.cells.volumes)./(trDNC), nlc * G.griddim);
   fprintf('----Compute alpha\n')
   alpha = sparse((1:sum(nlc)*G.griddim)', rldecode((1:G.cells.num)', nlc*G.griddim), 1) * ...
           ((trD .* G.cells.volumes)./(trDNC));

   alpha = alpha .* alpha_scaling;
   
   fprintf('----Compute SE\n')
   SE_AD = SparseMultiArray(alpha, [(1:nldofs)', (1:nldofs)'], {'i', 'j'}); 
   SE = SE_AD.asMatrix({'i', 'j'});

   fprintf('----Compute D\n')
   D_all = D.asMatrix({{'cell'}, {'i', 'j'}});
   D_AD = D;
   D = reshape(D_all(cells, :)', nlin, [])';
   [i, j] = blockDiagIndex(nlin * ones(numel(cells), 1), nlin * ones(numel(cells), 1));
   fprintf('----Final step of computing D\n');
   D = sparse(i, j, reshape(D', [], 1), nlin * numel(cells), nlin * numel(cells));
   fprintf('----computing KH\n');
   KH = op.volmap * op.WC * D * op.WC' + (op.I - op.PP)' * SE * (op.I - op.PP);
   
   fprintf('----computing S\n');
   S = op.assemb * KH * op.assemb';
    
   % adding computed matrices to return structure
   fprintf('----preparing return structure.\n');
   op.trD = trD;
   op.SE = SE;
   op.SE_AD = SE_AD;
   op.D = D;
   op.D_AD = D_AD;
   op.KH = KH;
   fprintf('---- all set\n');
end
