function [G, order, f] = glue2DGrid(G1, G2)
%Connect two 2D grids along common edge
%
% SYNOPSIS:
%    G        = glue2DGrid(G1, G2)
%   [G, iMap] = glue2DGrid(G1, G2)
%
% PARAMETERS:
%   G1    - First grid to be combined.
%
%   G2    - Second grid to be combined.
%
% DESCRIPTION:
%           Grids must follow definition from 'grid_structure'.  Both input
%           grids must be strictly two-dimensional both in terms of
%           `griddim` and in terms of size(nodes.coords, 2).
%
% RETURNS:
%   G    - Resulting grid structure.  Does not contain any Cartesian
%          information.  In particular, neither the `cells.indexMap` nor
%          the `cartDims` fields are returned--even when present in both
%          input grids.
%
%          Empty array ([]) if the input grids do not have a common
%          (non-empty) intersecting edge.
%
%   iMap - Input order map.  Two-element vector containing strictly either
%          [1,2] or [2,1].  This is the order in which the input grids (G1
%          and G2) are connected along the common intersecting edge.  If
%          this is [1,2] then the grids are geometrically concatenated as
%          [G1, G2].  Otherwise, the grids are geometrically concatenated
%          as [G2, G1].
%
%          One possible application of `iMap` is to determine the order in
%          which to extract the input grid` `indexMap` arrays for purpose
%          of creating the `indexMap` property for the combined grid `G`.
%
%          Empty array ([]) if the input grids do not have a common
%          (non-empty) intersecting edge.
%
% NOTE:
%   The result grid (`G`) does not provide derived geometric primitives
%   (e.g., cell volumes).  Such information must be explicitly computed
%   through a subsequent call to function `computeGeometry`.
%
%   To that end, the final step of function glue2DGrid is to order the
%   faces of all cells in a counter-clockwise cycle.  This sorting process
%   guarantees that face areas will be positive and that all face normals
%   will point from the first to the second cell in `G.faces.neighbors`.
%
%   This sorting step however is typically quite expensive and there are,
%   consequently, practical restrictions on the size of the input grids
%   (i.e., in terms of the number of cells) to function `glue2DGrid`.
%
% SEE ALSO:
%   `computeGeometry`

%{
Copyright 2009-2018 SINTEF ICT, Applied Mathematics.

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

   check_input(G1, G2);

   [bb{1:2}]  = deal(bounding_box(G1), bounding_box(G2));
   [lft, abo] = deal(left(bb{:}), above(bb{:}));

   if xor(lft ~= 0, abo ~= 0),
      % G1 either left/right or above/below G2.  Proceed accordingly.
      assert (sum(sum(bsxfun(@eq, [lft; abo], [1, 2]))) == 1, ...
              'Internal error.');

      if lft ~= 0,
         [connect, ord] = deal(@glue_left , lft);
      else
         [connect, ord] = deal(@glue_above, abo);
      end

      t     = { G1 , G2 };
      order = [ ord, 3 - ord ];  % 3 - ord == (2 - ord) + 1
      [G,f]     = connect(t{order});

      G.type = { mfilename };
      G_old=G;
      % get read of extra nodes HACK???
      [coords,ii,jj]=unique(G_old.nodes.coords,'rows');
      G=G_old;
      G.nodes.coords=coords;
      G.nodes.num=numel(ii);
      G.faces.nodes=jj(G.faces.nodes);
   else
      % G1 is neither left/right nor above/below G2 or BOTH of those
      % conditions simultaneously satisfied (single common point?).
      %
      % This means that no geometric connection between the two grids is
      % possible.  Return an empty result to indicate glue failure.
      G     = [];
      order = [];

      if mrstVerbose,
         n1 = inputname(1);  if isempty(n1), n1 = 'G1'; end
         n2 = inputname(2);  if isempty(n2), n2 = 'G2'; end

         fprintf(['Grids ''%s'' and ''%s'' don''t have a common', ...
                  'non-empty edge. Connection impossible.\n'], n1, n2);
      end
   end
end

%--------------------------------------------------------------------------

function check_input(G1, G2)
   assert (all([G1.griddim, G2.griddim] == 2) && ...
           all([size(G1.nodes.coords, 2), ...
                size(G2.nodes.coords, 2)] == 2) && ...
           all(diff(G1.faces.nodePos) == 2) && ...
           all(diff(G2.faces.nodePos) == 2), ...
           'Function ''%s'' is only supported in two space dimensions', ...
           mfilename);

   if (size(G1.cells.faces, 2) == 1) || (size(G2.cells.faces, 2) == 1),
      error(msgid('cfDirection:Missing'), ...
           ['Function ''%s'' requires direction information for ', ...
            'cell faces.'], mfilename);
   end
end

%--------------------------------------------------------------------------

function bb = bounding_box(G)
   c  = G.nodes.coords;
   bb = [ min(c, [], 1) ; max(c, [], 1) ];
end

%--------------------------------------------------------------------------

function c = left(bb1, bb2)
   c = position_impl(bb1, bb2, 1);
end

%--------------------------------------------------------------------------

function c = above(bb1, bb2)
   c = position_impl(bb1, bb2, 2);
end

%--------------------------------------------------------------------------

function c = position_impl(bb1, bb2, col)
   m = [ bb1(1, col) , bb2(1, col) ];  % Minimum coordinate
   M = [ bb1(2, col) , bb2(2, col) ];  % Maximum coordinate

   x = [ M(1) == m(2) , ...  % bb1 "left" of bb2
         m(1) == M(2) ];     % bb2 "left" of bb1

   c = find(x);

   if numel(c) ~= 1,
      % Neither is "left" of the other or both are "left" of the other.
      c = 0;
   end
end

%--------------------------------------------------------------------------

function [G,f] = glue_left(G1, G2)
% G1 left of G2.  Glue along common vertical edge.
   tag = [ 2, 1 ]; % Type 2 on left, 1 on right
   col = 2;        % Vertical edge => column 2

   [G,f] = glue_impl(G1, G2, tag, col);
end

%--------------------------------------------------------------------------

function [G,f] = glue_above(G1, G2)
% G1 above G2.  Glue along common horizontal edge
   tag = [ 4, 3 ]; % Type 4 on upper, 3 on lower
   col = 1;        % Horizontal edge => column 1

   [G,f] = glue_impl(G1, G2, tag, col);
end

%--------------------------------------------------------------------------


function [G,f] = glue_impl(G1, G2, tag, col)
   [f1, n1, i1, x1] = select_bfaces(G1, tag(1));
   [f2, n2, i2, x2] = select_bfaces(G2, tag(2));

   common   = intersection(x1, x2, col);
   affected = @(x) between(common(1), common(2), x(:, col));

   fnod = @(G, f)       G.faces.nodes(node_indices(G, f));
   lfe  = @(G, f, i)    i(reshape(fnod(G, f), 2, []) .');
   elim = @(G, f, i, a) f(sum(a(lfe(G, f, i)), 2) > 0);

   % Identify all nodes along common intersection.
   fe1 = elim(G1, f1, i1, affected(x1));
   il  = i1(fnod(G1, fe1));  xl = x1(il, col);

   fe2 = elim(G2, f2, i2, affected(x2));
   ir  = i2(fnod(G2, fe2));  xr = x2(ir, col);

   [N, ii] = intersection_topology(G1, G2, fe1, fe2, xl, xr);
   inodes  = intersection_geometry(G1.nodes.num, ii, n1, n2, il, ir);

   % Remove existing grid faces affected by common intersection.
   H1 = removeFaces(G1, fe1);
   H2 = removeFaces(G2, fe2);

   G = concat_grids(H1, H2);

   if col==1
      in = flipud(reshape(inodes,2,[])); inodes = in(:);
   end

   [G,f] = addFaces(G, inodes, repmat(2, [size(N, 1), 1]), N);
end

%--------------------------------------------------------------------------

function [N, ii] = intersection_topology(G1, G2, f1, f2, xl, xr)
   mby2 = @(a) reshape(a, 2, []) .';

   [u, ii, iu] = unique([xl ; xr], 'first');

   left     = mby2(iu(1 : numel(xl)));
   [ol, ol] = sort(left(:,1));                                  %#ok<ASGLU>
   left     = left(ol, :);

   right    = mby2(iu((numel(xl) + 1) : end));
   [or, or] = sort(right(:,1));                                 %#ok<ASGLU>
   right    = right(or, :);

   merge = NaN([numel(u) - 1, 2]);
   pl    = 1;
   pr    = 1;

   ix = diff(right,1,2)<0; right(ix,:) = right(ix,[2 1]);
   ix = diff(left, 1,2)<0; left (ix,:) = left (ix,[2 1]);

   for ival = 1 : (numel(u) - 1),
      % Left pointer
      [merge, pl] = record_interval(merge, pl, ival, 1, ol, left(pl,:));

      % Right pointer
      [merge, pr] = record_interval(merge, pr, ival, 2, or, right(pr,:));
   end

   bf2cell = @(G, f) sum(double(G.faces.neighbors(f, :)), 2);

   N = merge;
   i = ~isnan(N(:,1));
   N(i,1) = bf2cell(G1, f1(N(i,1)));

   i = ~isnan(N(:,2));
   N(i,2) = G1.cells.num + bf2cell(G2, f2(N(i,2)));

   N(isnan(N)) = 0;
end

%--------------------------------------------------------------------------

function inodes = intersection_geometry(nn1, ii, n1, n2, i1, i2)
   n1 = find(n1);
   n2 = find(n2);

   ni1 = numel(i1);
   p1  = ii <= ni1;

   inodes = zeros(size(ii));
   inodes(  p1) = n1(i1(ii(  p1)));
   inodes(~ p1) = n2(i2(ii(~ p1) - ni1)) + nn1;

   % Double internal nodes to account for beginning *and* end of intervals.
   inodes = inodes(1 + fix((1 : 2*(numel(inodes) - 1)) ./ 2));
end

%--------------------------------------------------------------------------

function [f, n, i, x] = select_bfaces(G, tag)
   bf        = boundary_faces(G);
   cf        = G.cells.faces(bf.i, :);
   f         = cf(bf.e(cf(:,1)) & (cf(:,2) == tag), 1);
   %f         = cf(bf.e(cf(:,1)), 1);
   [n, i, x] = face_nodes(G, f);
end

%--------------------------------------------------------------------------

function bf = boundary_faces(G)
   e = any(G.faces.neighbors == 0, 2);

   c = false([G.cells.num, 1]);
   c(sum(G.faces.neighbors(e, :), 2)) = true;
   c = find(c);

   i = mcolon(G.cells.facePos(  c  ), ...
              G.cells.facePos(c + 1) - 1);

   bf = struct('e', e, 'c', c, 'i', reshape(i, [], 1));
end

%--------------------------------------------------------------------------

function [n, i, x] = face_nodes(G, f)
   n = false([G.nodes.num, 1]);
   i = NaN  (size(n));

   p = node_indices(G, f);

   n(G.faces.nodes(p)) = true;

   i(n) = 1 : sum(n);
   x    = G.nodes.coords(n, :);
end

%--------------------------------------------------------------------------

function i = node_indices(G, f)
   i = mcolon(G.faces.nodePos(f), G.faces.nodePos(f + 1) - 1);
end

%--------------------------------------------------------------------------

function common = intersection(x1, x2, col)
   lo = @(x) min(x(:, col), [], 1);
   hi = @(x) max(x(:, col), [], 1);

   common = [ max(lo(x1), lo(x2)) , ...
              min(hi(x1), hi(x2)) ];

   assert (common(1) < common(2), 'Non-empty intersection?');
end

%--------------------------------------------------------------------------

function tf = between(lo, hi, x)
   assert (lo < hi, 'Internal error');
   tf = ~(x < lo) & ~(hi < x);
end

%--------------------------------------------------------------------------

function [merge, k] = record_interval(merge, k, ival, col, o, x)
   [lo, hi] = deal(x(1), x(2));

   if (lo <= ival) && (ival < hi),
      % ival subset of [lo, hi).  Record that fact.
      merge(ival, col) = o(k);
   end

   if (hi == ival + 1) && (k < numel(o)),
      k = k + 1;
   end
end

%--------------------------------------------------------------------------

function G = concat_grids(G1, G2)
   cells = concat_cells(G1, G2);
   faces = concat_faces(G1, G2);
   nodes = concat_nodes(G1, G2);

   G = struct('nodes'  , nodes, ...
              'cells'  , cells, ...
              'faces'  , faces, ...
              'griddim', 2);
end

%--------------------------------------------------------------------------

function cells = concat_cells(G1, G2)
   fpos    = concat_pos(G1.cells.facePos, G2.cells.facePos);
   fadd    = zeros([1, size(G2.cells.faces, 2)]);
   fadd(1) = G1.faces.num;
   faces   = [G1.cells.faces; bsxfun(@plus, G2.cells.faces, fadd)];

   num  = G1.cells.num + G2.cells.num;
   imap = []; % Concatenated grids have no sensible indexMap.

   %tag=[];
   if(isfield(G1.cells,'tag'))
       tag=G1.cells.tag;
   else
       tag=zeros(G1.cells.num,1);
   end
   if(isfield(G2.cells,'tag'))
       tag=[tag;G2.cells.tag];
   else
       tag=[tag;zeros(G2.cells.num,1)];
   end
   cells = struct('num'     , num  , ...
                  'facePos' , fpos , ...
                  'faces'   , faces, ...
                  'indexMap', imap,...
                  'tag',tag);
end

%--------------------------------------------------------------------------

function faces = concat_faces(G1, G2)
   num   = G1.faces.num + G2.faces.num;
   npos  = concat_pos(G1.faces.nodePos, G2.faces.nodePos);
   nods  = [G1.faces.nodes; G1.nodes.num + G2.faces.nodes];

   remap = [ 0 ; G1.cells.num + (1 : G2.cells.num).' ];
   neigh = [double(G1.faces.neighbors); ...
            remap(G2.faces.neighbors + 1)];

   faces = struct('num'      , num  , ...
                  'nodePos'  , npos , ...
                  'nodes'    , nods , ...
                  'neighbors', neigh);

   tag = [];
   if isfield(G1.faces, 'tag') || isfield(G2.faces, 'tag'),
      tag = zeros([faces.num, 1]);

      if isfield(G1.faces, 'tag'),
         tag(1 : G1.faces.num)         = G1.faces.tag;
      end

      if isfield(G2.faces, 'tag'),
         tag((G1.faces.num + 1) : end) = G2.faces.tag;
      end
   end

   if ~isempty(tag),
      faces.tag = tag;
   end
end

%--------------------------------------------------------------------------

function nodes = concat_nodes(G1, G2)
   nodes = struct('num'   , G1.nodes.num + G2.nodes.num, ...
                  'coords', [G1.nodes.coords; G2.nodes.coords]);
end

%--------------------------------------------------------------------------

function pos = concat_pos(p1, p2)
   pos = [ p1(1 : (end - 1)) ; cumsum([ p1(end) ; diff(p2) ]) ];
end
