function [V, P, Q] = evalWellBasis(well, g, cg, BI, C, D, w, mob, varargin)
%Compute multiscale basis functions for single well.
%
% SYNOPSIS:
%   [V, P, Q] = evalWellBasis(w, G, CG, BI, C, D, weight, mob)
%   [V, P, Q] = evalWellBasis(w, G, CG, BI, C, D, weight, mob, ...
%                             'pn1', pv1, ...)
%
% PARAMETERS:
%   w        - Single well structure as defined by functions 'addWell' and
%              'assembleWellSystem'.
%
%   G, CG    - Underlying grid (G, see grid_structure) and coarse grid
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
%   'pn'/pv  - List of 'key'/value pairs defining optional parameters.  The
%              supported options are:
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
%                - OverlapWell  --
%                         Number of fine grid cells by which to extend
%                         support of the well basis function about the well
%                         bore.
%
%                - OverlapBlock --
%                         Number of fine grid cells by which to extend
%                         support of the well basis function about the
%                         reservoir coarse block.
%
%                - LinSolve --
%                         Handle to linear system solver software to which
%                         the fully assembled system of linear equations
%                         will be passed.  Assumed to support the syntax
%
%                             x = LinSolve(A, b)
%
%                         in order to solve a system Ax=b of linear eq.
%                         Default value: LinSolve = @mldivide (backslash).
%
% RETURNS:
%   V - Cell array, one element for each coarse block intersected by the
%       well 'w', of cell arrays of SPARSE input vectors from which the
%       flux basis function matrix B*\Psi may be formed.
%
%   P - Cell array, one element for each coarse block intersected by the
%       well 'w', of cell arrays of SPARSE input vectors from which the
%       pressure basis function matrix Phi may be formed.
%
%   Q - Matrix of size nwc-by-ncb of (well-to-block) well rate basis
%       functions.  Here, 'nwc' is the number of perforations (i.e.,
%       cells intersected by) the well 'w'.
%
% SEE ALSO:
%   `evalBasisFunc`, `generateCoarseWellSystem`.

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


   opt = struct('src', [], 'OverlapWell', 0, 'OverlapBlock', 0, ...
                'LinSolve',     @mldivide);
   opt = merge_options(opt, varargin{:});

   part    = get_partition(cg);
   sub_map = mappings(g, part, well, mob, opt);
   theta   = get_coarse_weighting(g, part, w, opt.src);
   avgmob  = average_mobility(part, mob, g.cells.volumes);

   pwc = part(well.cells);
   blk = unique(pwc);

   iq     = []; jq = []; sq = [];
   [V, P] = deal(cell([numel(blk), 1]));

   nBlk     = max(part);
   src_mult = zeros([nBlk       , 1]);
   press_ok = zeros([numel(blk) , 1]);
   renum    = zeros([g.cells.num, 1]);

   pick = false([g.cells.num, 1]);
   pick(well.cells) = true;

   for k = 1 : numel(blk),
      b  = blk(k);
      iG = sub_map.cells(b) ;   % Fine-scale cells      present in 'b'

      % hack to get rid of well cells that are not connected to the block
      c = min(iG(~pick(iG)));
      if isempty(c),
         iG=iG(part(iG)==b);
      else
         p = ones(g.cells.num,1); p(iG)=2;
         p = processPartition(g,p);
         iG = find(p==p(c));
      end

      iF = sub_map.hf   (iG);   % Fine-scale half-faces present in 'b'
      iH = sub_map.faces(iF);   % Fine-scale faces      present in 'b'

      renum(iG) = 1 : numel(iG);

      % Extract linsys matrix components for sub-(hf,cells,faces) in
      % reservoir part of model.
      sBI  = BI(iF,iF); sC = C(iF,iG); sD = D(iF,iH);

      % Build reservoir linsys right hand side components for hybrid system.
      %  -1 source strength in b
      %   0 source strength elsewhere (if opt.Overlap{Well,Block} > 0)
      %
      src_mult(b) = -1;
      sG          = theta(iG) .* src_mult(part(iG));
      src_mult(b) = 0;

      sF = zeros(size(iF)); sH = zeros(size(iH));

      % Append well block 'b' well system to reservoir discrete system.
      iFw  = find(pwc == b);
      niFw = numel(iFw);
      wc   = well.cells(iFw);
      dmob = blkdiag(sub_map.dmob(iF,numel(iF)),       ...
                     spdiags(mob(wc), 0, niFw, niFw));

      sBI = dmob * blkdiag(sBI, spdiags(well.WI(iFw), 0, niFw, niFw));
      sC  = vertcat(sC, sparse(1 : niFw, renum(wc), 1, niFw, size(sC,2)));
      sD  = blkdiag(sD, sparse(niFw, 0));
      sF  = vertcat(sF, zeros([niFw, 1]));  % Model as bhp w/target p=0.

      % Solve composite system of linear equations.
      [v, p, lam] = schurComplementSymm(sBI, sC, sD, sF, sG, sH, ...
                                        'LinSolve', opt.LinSolve);

      Bv = dmob * ([sC, sD] * [p; -lam]);  % Bv <- B*v == C*p - D*lam

      % Orthogonalize pressure against source term (theta(iG)) in each
      % coarse block.
      [p, press_ok(k)] = orth_press(p, theta(iG), part(iG), nBlk);

      j1 = zeros(size(iF ));  j1(1) = 1;
      j2 = zeros(size(iG ));  j2(1) = 1;
      j3 = zeros(size(iFw));  j3(1) = 1;

      % SPARSE input triplets [i,j,v] used in unpackWellSystemComponentsMS,
      % by way of solveIncompFlowMS>syscomp_wells, to form flux- and
      % pressure basis function matrices [Psi and Phi, respectively].  The
      % average mobility is, strictly speaking, only needed to form 'Phi',
      % but we store the value for 'Psi' as well.  If memory becomes a
      % (serious) issue, we may decide otherwise.
      %
      V{k} = {iF, j1, Bv(1:numel(iF)), avgmob(b)};
      P{k} = {iG, j2, p              , avgmob(b)};

      iq = [iq; iFw]; jq = [jq; j3];      %#ok
      sq = [sq; -v(numel(iF) + 1 : end)]; %#ok
   end

   jq = cumsum(jq);
   Q  = sparse(iq, jq, sq, numel(pwc), jq(end));

   if ~all(press_ok),
      warning(msgid('Pressure:Orthogonality'), ...
             ['At least one pressure basis function does not ', ...
              'satisfy orthogonality condition.\nQuality of ' , ...
              'solution may be reduced.']);
   end
end

%--------------------------------------------------------------------------
% Helpers follow.
%--------------------------------------------------------------------------

function p = get_partition(cg)
   p = cg.partition;
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

function ops = mappings(g, p, well, mob, opt)
   cellNo = rldecode(1 : g.cells.num, diff(g.cells.facePos), 2) .';

   nc = g.cells.num;
   N  = double(g.faces.neighbors(all(g.faces.neighbors > 0, 2), :));
   N  = sparse([N(:,1); N(:,2); (1 : nc).'], ...
               [N(:,2); N(:,1); (1 : nc).'], 1, nc, nc);

   wc = sparse(well.cells, 1, 1, nc, 1);
   for k = 1 : opt.OverlapWell, wc = N * wc; end
   wc = logical(wc);

   %-----------------------------------------------------------------------
   % Define mapping operators. --------------------------------------------
   %
   %   1) ops.cells:
   %      any(sub_c(:,b), 2) is true for all fine-scale cells in the
   %      (possibly extended) block(s) 'b'.
   %
   ops.cells = @(b) find(any([sub_c(N, p==b, opt.OverlapBlock), wc], 2));

   %   2) ops.hf:
   %      Transpose because MCOLON gives row-vector, while we need columns.
   %
   ops.hf    = @(c) mcolon(g.cells.facePos(  c  ), ...
                           g.cells.facePos(c + 1) - 1) .';

   %   3) ops.faces:
   %      FIND those faces mentioned at least once in the half-faces 'i'.
   %
   ops.faces = @(i) find(accumarray(g.cells.faces(i,1), 1) > 0);

   %   4) ops.dmob:
   %      mob(cellno(i)) is (total) mobility in cell containing half-face
   %      'i'.
   %
   ops.dmob  = @(i,n) spdiags(mob(cellNo(i)), 0, n, n);
end

%--------------------------------------------------------------------------

function c = sub_c(N, c0, overlap)
   c     = zeros(size(c0));
   c(c0) = 1;

   for k = 1 : overlap, c = N * c; end

   c = sparse(logical(c));
end

%--------------------------------------------------------------------------

function mob = average_mobility(p, mob, vol)
   assert (all(numel(p) == [numel(mob), numel(vol)]));

   mob = accumarray(p, mob .* vol) ./ accumarray(p, vol);
end

%--------------------------------------------------------------------------

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
