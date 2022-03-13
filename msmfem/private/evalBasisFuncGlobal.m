function [V, P] = evalBasisFuncGlobal(faces, g, cg, BI, C, D, w, mob, ...
                                      global_inf, varargin)
%Compute multiscale basis functions for selected faces.
%
% SYNOPSIS:
%   [V, P] = evalBasisFunc(faces, G, CG, BI, C, D, weight, mob)
%   [V, P] = evalBasisFunc(faces, G, CG, BI, C, D, weight, mob, ...
%                          'pn1', pv1, ...)
%
% PARAMETERS:
%   faces    - List (array) of coarse faces for which (new) basis function
%              values are requested.  Must be a list of explicit indices.
%              Function 'evalBasisFunc' does not support a LOGICAL array.
%              I: NB: how to choose which faces to add? now that we might
%              have global info for all faces.. even 0-neumann.. should
%              0-neumann faces be exluded?! probably!!
%
%   G, CG    - Underlying grid (G, see 'grid_structure') and coarse grid
%              model (CG) defined by function 'generateCoarseGrid'.
%
%   BI, C, D - Component system matrices from underlying fine-grid
%              (mimetic) hybrid discretisation of incompressible pressure
%              equation.  Note that 'BI' is assumed to be the hybrid system
%              INVERSE mass matrix (and SPD).
%
%   weight   - Array of numerically evaluated synthetic volume source
%              weighting term.  One scalar value for each cell in the
%              underlying fine grid G.  All values must be supplied, even
%              if some of the cells do not participate in the any of the
%              basis functions for 'faces'.
%
%   mob      - Total mobility.  One scalar value for each cell in the
%              underlying fine model.
%
%   global_inf - global information from fine scale solution
%                (fineSol.faceFlux) to be used as boundary condition when
%                calculating basis for coarse faces in the interior of the
%                domain.
%
%   'pn'/pv  - List of 'key'/value pairs defining optional parameters.  The
%              supported options are:
%                - Verbose --
%                         Whether or not to display and update a progress
%                         bar whilst computing the basis functions.
%                         Logical.  Default value depending on the global
%                         verbose settings of function 'mrstVerbose'.
%
%                - src -- Explicit source terms in the underlying fine grid
%                         model which must be taken into account when
%                         generating the basis functions.  Note that basis
%                         functions in blocks containing an explicit source
%                         term will be generated based solely on the
%                         explicit source.  The values in 'weight'
%                         pertaining to such blocks will be ignored.
%
%                         Must be a source data structure as defined by
%                         function 'addSource'.  Default value is [] (an
%                         empty array), meaning that no explicit sources
%                         are present in the model.
%
%                - bc  -- External boundary conditions.
%                         Must be a boundary condition data structure as
%                         defined by function 'addBC'.  Default value is []
%                         (an empty array), meaning that no boundary
%                         conditions are present in the model.
%
%                - Overlap --
%                         The number of fine-scale cells by which to extend
%                         the support of a given flux or pressure basis
%                         function into neighbouring coarse blocks.  Note
%                         well: Using Overlap > 0 precludes hybrid
%                         formulation of the resulting coarse system of
%                         linear equations.  Default value: Overlap = 0
%                         (don't extend the basis function support).
%
%                - LinSolve --
%                         Handle to linear system solver software to which
%                         the fully assembled system of linear equations
%                         will be passed.  Assumed to support the syntax
%
%                         x = LinSolve(A, b)
%
%                         in order to solve a system Ax=b of linear eq.
%                         Default value: LinSolve = @mldivide (backslash).
%
% RETURNS:
%   V - Cell array, one element for each coarse face in 'faces', of cell
%       arrays (tuples) of SPARSE input vectors (and auxillary information)
%       from which the flux basis function matrix B*\Psi may be formed.
%
%   P - Cell array, one element for each coarse face in 'faces', of cell
%       arrays (tuples) of SPARSE input vectors (and auxillary information)
%       from which the pressure basis function matrix Phi may be formed.
%
% EXAMPLE:
%   % Define simple geometry and rock properties.
%   G = computeGeometry(cartGrid([50, 50, 4]));
%   rock = struct('perm' , repmat(0.1*darcy(), [G.cells.num, 1]));
%
%   % Define coarse grid.
%   CG = generateCoarseGrid(G, partitionUI(G, [5, 5, 1]));
%
%   % Assemble hybrid pressure system matrices.
%   S = computeMimeticIP(G, rock);
%
%   % Compute basis functions for all inner coarse faces using cell
%   % permeability as a weighting function.
%   %
%   faces  = find(all(CG.faces.neighbors > 0, 2));
%   [V, P] = evalBasisFunc(faces, G, CG, S.BI, S.C, S.D, ...
%                          rock.perm, ones([G.cells.num, 1]));
%
% SEE ALSO:
%   `generateCoaseSystem`, `generateCoarseGrid`, `assignBasisFuncs`,
%   `computeMimeticIP`, `mrstVerbose`.

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


   assert( size(global_inf,1) == g.faces.num);
   if ~isnumeric(faces),
      error(msgid('Faces:NonNumeric'), ...
            'Argument ''faces'' must be a numeric array.');
   end

   opt = struct('Verbose', mrstVerbose, 'Overlap', 0, ...
                'src', [], 'bc', [], ...
                'LinSolve',   @mldivide);

   opt = merge_options(opt, varargin{:});

   part = get_partition(cg);
   [sub_map, ncells, nhf_blk] = mappings(g, cg, part, mob, global_inf, opt);
   nBlk = max(part);

   avgmob = average_mobility(part, mob, g.cells.volumes);

   % Driving forces (sources and BCs).
   [f_bc, h_bc, is_dirichlet] = expand_bc(g, opt.bc);
   theta = get_coarse_weighting(g, part, w, opt.src);

   if opt.Verbose,
      h = waitbar(0, ...
                 ['Computing flux and pressure basis functions ', ...
                  'using global information...']);
      nFaces = numel(faces);
      fprintf(['Computing flux and pressure basis functions ', ...
               'using global information... ']), tic
   end

   V = cell([numel(faces), 1]);
   P = cell([numel(faces), 1]);

   sgn      = [1, -1];
   press_ok = false([numel(faces), 1]);
   for f = 1 : numel(faces),
      face = faces(f);
      blk  = cg.faces.neighbors(face,:); % Blocks connected to 'face'.
      blk  = blk(blk > 0);

      nblk = numel(blk);                 % Number of blocks sharing 'face'.

      v_  = cell([nblk, 1]);
      p_  = cell([nblk, 1]);
      iF_ = cell([nblk, 1]);
      iG_ = cell([nblk, 1]);

      % The basis for the coarse face 'face' is a combination of the bases
      % in blk. For each each block in blk compute basis by using boundary
      % conditions from either
      % 1) the fine system - if 'face' is on the boundary,
      % 2) the fine-scale solution - if 'face' is in the int. of the domain.
      for i = 1:nblk
         block = blk(i);

      iG = sub_map.cells(block); % Fine-scale cells      present in 'block'
      iF = sub_map.hf   (iG) ;   % Fine-scale half-faces present in 'block'
      iH = sub_map.faces(iF) ;   % Fine-scale faces      present in 'block'

      % Extract linsys matrix components for sub-(hf,cells,faces).
      sBI = BI(iF,iF); sC = C(iF,iG); sD = D(iF,iH);

      % Build linsys right hand side components for hybrid system.
      % 1) Define driving source term:
      %      +1 source strength in blk(1)
      %      -1 source strength in blk(2) (if nblk > 1)
      %       0 source strength elsewhere (if opt.Overlap > 0)
      %
      src_mult = sgn(i);
      sG       = theta(iG).* src_mult;

      % 2) Define trivial vectors [f,h] (no external forces).
      %    Will be non-trivial in case of boundary conditions (below).
      niF = numel(iF);
      sF  = zeros([niF, 1]); sH = zeros([numel(iH), 1]);

      % Update linsys components for presence of several phases...
      dmob = sub_map.dmob(iF, niF);
      sBI  = dmob * sBI;

      is_dir = false(size(sH));
      lam    = zeros(size(sH));

      if nblk == 1,
         % The coarse face 'face' is on the boundary of the domain.
         % Consequently, we need to include boundary conditions for 'face'
         % but not for any other external coarse faces which might happen
         % to connect to 'blk'.  Such other faces always get no-flow
         % conditions (homogeneous Neumann).
         %
         [sF, sH, lam, is_dir] = handle_bc(sF, sH, lam, is_dir,   ...
                                           g, f_bc, h_bc, iF, iH, ...
                                           is_dirichlet,          ...
                                           sub_map.sub_f(face));

      else
         % The coarse face 'face' is inside the domain and we use the
         % global info as boundary condition to compute the blk(i)-part of
         % the basis for face.
         sH = handle_global(sH, g, global_inf, iH, ...
                                    sub_map.sub_f(face), block, part, i);
      end

      do_reg = ~any(is_dir);  % Need to set pressure zero level?
      [v, p, lam(~is_dir)] = schurComplementSymm(sBI, sC, sD(:,~is_dir), ...
                                                 sF , sG, sH(  ~is_dir), ...
                                                 'Regularize', do_reg, ...
                                                 'LinSolve', opt.LinSolve);

      v_{i} = dmob * ([sC, sD] * [p; -lam]);   % v <- B*v == C*p - D*lam

      % Orthogonalize pressure against source term (theta(iG)) in each
      % coarse block.
      [p_{i}, press_ok(f)] = orth_press(p, theta(iG), part(iG), nBlk);
      iF_{i} = iF;
      iG_{i} = iG;
      end

      % entries corresponding to blk(1) stored *first* in the following
      % SPARSE tuples.

      % SPARSE input duplets [i,s] and auxillary information used in
      % solveIncompFlowMS>syscomp_res to form either the hybrid or mixed
      % flux and pressure basis function matrices 'Psi' and 'Phi',
      % respectively.  The tuples are
      %
      %    1  2  3  4  5  6  7
      %   {i, s, f, b, n, m, o}
      %
      % with 'f' being the coarse face to which the values are associated,
      % 'b' being the coarse block(s) connected to 'f'.  Furthermore, 'n'
      % is the number of entries in [i,s] belonging to the *first* coarse
      % block (b(1)) while 'm' and 'o' are the average mobility in block(s)
      % 'b' and number of fine-scale cells by which the support of a basis
      % function is extended about the primary domain, 'b'.  The average
      % mobility is used in the formation of 'Phi'
      % (solveIncompFlowMS>syscomp_res) while the amount of overlap is used
      % only to determine the suitability of the coarse hybrid solver. The
      % coarse hybrid solver cannot be employed if any basis function is
      % generated with a positive amount of overlap.
      %
      % The entry 'n' (i.e., V{k}{5} and P{k}{5}) is only used when forming
      % the *hybrid* basis matrices (solveIncompFlowMS>split_for_hybrid).
      % We note that the remaining SPARSE input vector 'j' is formed
      % differently depending on which type of basis function matrix to
      % construct (hybrid or mixed).  See the sub functions
      % split_for_hybrid and expand_mixed in solveIncompFlowMS for details.
      %

      if opt.Verbose, waitbar(f / nFaces, h), end

      V{f} = {vertcat(iF_{:}), vertcat(v_{:}), face, blk, ...
              nhf_blk(blk), avgmob(blk), opt.Overlap};
      P{f} = {vertcat(iG_{:}), vertcat(p_{:}), face, blk, ...
              ncells(blk) , avgmob(blk), opt.Overlap};
   end
   if opt.Verbose,
      toc, close(h)
      if ~all(press_ok),
         warning(msgid('Pressure:Orthogonality'), ...
                ['At least one pressure basis function does not ', ...
                 'satisfy orthogonality condition.\nQuality of ' , ...
                 'solution may be reduced.']);
      end
   end
end


%--------------------------------------------------------------------------
% Helpers follow.
%--------------------------------------------------------------------------


function p = get_partition(cg)
% get_partition -- Extract explicit partition vector from coarse grid.
%
% SYNOPSIS:
%   p = get_partition(CG)
%
% PARAMETERS:
%   CG - Coarse grid structure.
%
% RETURNS:
%   p  - Partition vector used in creating the coarse grid.

   p = cg.partition;
end

%--------------------------------------------------------------------------

function [ops, nc, nhf] = mappings(g, cg, p, mob, global_inf, opt)
% mappings - Build essential grid mappings for coarse grid.
%
% SYNOPSIS:
%   [ops, nc, nhf] = mappings(G, CG, p, mob, opt)
%
% PARAMETERS:
%   G, CG - Fine grid and coarse grid data structures, respectively.
%   p     - Partition vector, size [G.cells.num,1], such that p(i) is the
%           coarse block containing cell 'i'.  This is assumed to be the
%           original partition vector which created the coarse grid.
%   mob   - Vector, size [G.cells.num,1], of total mobility values for all
%           cells in the fine-scale model.
%   opt   - Option structure describing external influences (and, in
%           particular, external boundary conditions).  Also, assumed to
%           contain a scalar field 'Overlap' denoting the number of cells
%           by which to extend support of the basis functions into
%           neighbouring coarse blocks.  Note: Overlap > 0 precludes hybrid
%           formulation of the resulting coarse system.
%
% RETURNS:
%   ops - A data structure whose fields are functions defining particular
%         aspects of mapping topological information from the fine grid to
%         the coarse grid.  Specifically, the fields are
%            - cells - Function @(b) -> fine cells attached to blocks 'b'.
%            - hf    - Function @(c) -> fine half-faces in fine cells 'c'.
%            - faces - Function @(i) -> fine faces corresponding to hf 'i'.
%            - sub_f - Function @(f) -> constituent fine faces of coarse
%                      face 'f'.  Defined non-trivially only when there are
%                      any external (coarse) faces supporting flow (i.e.,
%                      not no-flow conditions) and for which (flux) basis
%                      functions must be computed.
%            - dmob  - Function @(i,n) -> n-by-n diagonal matrix of (total)
%                      mobility values for all half-faces 'i'.
%   nc  - Array of size [CG.cells.num, 1] containing the number of
%         fine-scale cells in each block.  Specifically, nc(b) is the
%         number of fine-scale cells in coarse block 'b'.
%   nhf - Array of size [CG.cells.num, 1] containing the number of
%         fine-scale half-faces in each block.  Specifically, nhf(b) is the
%         number of fine-scale half-faces in coarse block 'b'.

   % Compute 'sub_f' mapping depending on existence of any outer coarse
   % faces supporting flow (i.e., BC is not no-flow).  See 'help subFaces'
   % for details on MCOLON expression.
   %I:
   if ~isempty(opt.bc) || ~isempty(global_inf), %NB: always true
      [nsub, sub] = subFaces(g, cg);
      sub_ix      = cumsum([0; nsub]);
      sub_f       = @(f) sub(mcolon(sub_ix(f) + 1, sub_ix(f + 1)));
   else
      % Dangerous.  Assumes we won't be called upon to build a basis
      % function for an outer coarse face (f) if there are no external
      % boundary conditions or global info...
      sub_f       = @(f) [];
   end

   % Count number of cells and number of half-faces in all coarse blocks.
   nc  = accumarray(p, 1);
   nhf = accumarray(p, double(diff(g.cells.facePos)));

   % sub_c(:,b) == true for all cells (rows) within (extended) block 'b'.
   sub_c = sub_cells(g, p, opt.Overlap);

   % cellno(i) == (fine-scale) cell which contains half-face 'i'.
   % hfix      -- Index into 'G.cells.faces'.  Specifically, data concerning
   %              fine-scale cell 'c' occupies rows [hfix(c)+1 : hfix(c+1)]
   %              of 'G.cells.faces'.
   %
   cellno = rldecode(1 : g.cells.num, diff(g.cells.facePos), 2) .';

   %-----------------------------------------------------------------------
   % Define mapping operators. --------------------------------------------
   %
   %   1) ops.cells:
   %      any(sub_c(:,b), 2) is true for all fine-scale cells in the
   %      (possibly extended) block(s) 'b'.
   %
   ops.cells = @(b) find(any(sub_c(:,b), 2));

   %   2) ops.hf:
   %      Transpose because MCOLON gives row-vector, while we need columns.
   %
   ops.hf    = @(c) mcolon(g.cells.facePos(  c  ), ...
                           g.cells.facePos(c + 1) - 1) .';

   %   3) ops.faces:
   %      FIND those faces mentioned at least once in the half-faces 'i'.
   %
   % will be sorted..
   ops.faces = @(i) find(accumarray(g.cells.faces(i,1), 1) > 0);

   %   4) ops.sub_f:
   %      See above.
   %
   ops.sub_f = sub_f;

   %   5) ops.dmob:
   %      mob(cellno(i)) is (total) mobility in cell containing half-face
   %      'i'.
   %
   ops.dmob  = @(i,n) spdiags(mob(cellno(i)), 0, n, n);
end

%--------------------------------------------------------------------------

function theta = get_coarse_weighting(g, p, w, src)
   % All synthetic weighting terms must be strictly non-negative.
   %
   assert (~any(w < 0));

   % Initially, assume there are no explicit/external sources.
   %
   theta = w .* g.cells.volumes;

   if ~isempty(src),
      % Update for explicit sources if there nevertheless are some...

      % Determine coarse blocks already containing external sources.
      %
      has_src = accumarray(p(src.cell), 1, [max(p), 1]) > 0;

      % Eliminate previous (synthetic) weighting and apply correct source.
      %
      theta(has_src(p)) = 0;
      theta(src.cell)   = src.rate;
   end

   % Note:
   %   We need to normalize the (synthetic or explicit) source term 'theta'
   %   such that \int_{B_i} theta d\Omega == 1 lest the basis functions be
   %   inconsistent.
   %
   denom = accumarray(p, theta);  assert (all(abs(denom) > 0));
   theta = theta ./ denom(p);
end

%--------------------------------------------------------------------------

function [f, h, d] = expand_bc(g, bc)
   f = zeros([g.faces.num, 1]);
   h = zeros([g.faces.num, 1]);
   d = false([g.faces.num, 1]);

   if ~isempty(bc),
      assert (all(accumarray(bc.face, 1, [g.faces.num, 1]) <= 1));

      is_dir = strcmp('pressure', bc.type);
      f(bc.face(is_dir)) = bc.value(is_dir);
      d(bc.face(is_dir)) = true;

      is_neu = strcmp('flux', bc.type);
      h(bc.face(is_neu)) = - bc.value(is_neu);
   end
end

%--------------------------------------------------------------------------

function sub_c = sub_cells(g, p, overlap)
   nc    = g.cells.num;
   sub_c = sparse(1 : nc, p, 1, nc, max(p));  % == cg.cells.subCells

   if overlap > 0,
      n = double(g.faces.neighbors(all(g.faces.neighbors > 0, 2), :));
      n = sparse([n(:,1); n(:,2); (1 : nc).'], ...
                 [n(:,2); n(:,1); (1 : nc).'], 1, nc, nc);

      % BFS to discover immediate neighbours in overlap region.
      for o = 1 : overlap, sub_c = n * sub_c; end
   end

   sub_c = logical(sub_c);
end

%--------------------------------------------------------------------------

function mob = average_mobility(p, mob, vol)
   assert (all(numel(p) == [numel(mob), numel(vol)]));

   mob = accumarray(p, mob .* vol) ./ accumarray(p, vol);
end

%--------------------------------------------------------------------------
%orth_press(p, theta(iG), part(iG), nBlk);
function [p, ok] = orth_press(p, w, b, nb)
   present        = false([nb, 1]);
   present(b)     = true;

   renum          = zeros([nb, 1]);
   renum(present) = 1 : sum(present);

   % Compute orthogonality adjustment constant.  One value for each coarse
   % block present in the support of 'p'.
   a = accumarray(renum(b), p .* w) ./ ...
       accumarray(renum(b),      w);

   p = p - a(renum(b));

   % Check orthogonality to within relative bounds.
   ok = abs(p' * w) < 2 * numel(p) * eps(norm(p,inf));
end

%--------------------------------------------------------------------------

function [sF, sH, lam, is_dir] = handle_bc(sF, sH, lam, is_dir, ...
                                           g, f_bc, h_bc, iF, iH, ...
                                           is_dirichlet, ih)

   fno     = zeros([g.faces.num, 1]);
   fno(iH) = 1 : numel(iH);

   is_bdry          = false(size(sH));
   is_bdry(fno(ih)) = true;
   is_dir(:)        = is_bdry & is_dirichlet(iH);
   is_neu           = is_bdry & ~is_dir;

   % This code assumes that a coarse face is either Dirichlet or
   % Neumann (not both).  That assumption must be revisited if the
   % following assertion fails.  We also fail if the face is neither
   % Dirichlet nor Neumann.
   %
   tot_flx = norm(h_bc(iH(is_neu)));
   assert (xor(sum(is_dir) > 0 && tot_flx < 1000 * eps(1), ...
               ~(tot_flx < 1000 * eps(1)) && sum(is_dir) == 0));

   % Dirichlet (pressure) boundary conditions.
   lam(is_dir) = f_bc(iH(is_dir));
   loc_iF      = is_dir(fno(g.cells.faces(iF,1)));
   sF(loc_iF)  = sF(loc_iF) - lam(is_dir);

   % Neumann (flux) boundary conditions.
   sH(is_neu)  = h_bc(iH(is_neu));
   denom = sum(sH(is_neu));
   if abs(denom) > sqrt(eps(denom)),
      sH(is_neu) = sH(is_neu) ./ denom;
   end
end

%--------------------------------------------------------------------------

function [sH] = handle_global(sH, g, global_flux, iH, ih, blk, p, i)
   % merge with handle_bc?
   fno     = zeros([g.faces.num, 1]);
   fno(iH) = 1 : numel(iH);

   is_neu          = false(size(sH));
   is_neu(fno(ih)) = true;

   p1 = [0; p];
   outblk = p1(g.faces.neighbors(ih,1) + 1);

   cellSgn = 2*(outblk == blk)-1;

   % global_flux is given as faceFlux: the direction of the flux is from
   % blk(1) to blk(2). In hybrid systems flux bc values are given as
   % cellFlux. Thus, for blk(2) we must multiply sH by -1 to get correct
   % flow i.e. force direction of flow from blk(1) to blk(2).
   blockSgn = 1;
   if i == 2, blockSgn = -1; end

   % Neumann (flux) boundary conditions.
   sH(is_neu)  = cellSgn.*global_flux(ih);
   denom = sum(sH(is_neu));
   if abs(denom) > sqrt(eps(denom)),
      sH(is_neu) = blockSgn.*sH(is_neu) ./ denom;
   end
end
