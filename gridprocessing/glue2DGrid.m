function G = glue2DGrid(G1, G2, varargin)
%Connect two 2D grids along common edges
%
% SYNOPSIS:
%    G = glue2DGrid(G1, G2)
%    G = glue2DGrid(G1, G2, 'pn1', pv1, ...)
%
% PARAMETERS:
%   G1 - First grid to be combined.
%
%   G2 - Second grid to be combined.
%
% KEYWORD ARGUMENTS:
%   tol - Define geometric tolerance (Euclidian distance, metres) to
%         determine coincidence between nodes, and between nodes and edges,
%         between the two grids.  Default is 1e-3.
%
% LIMITATIONS:
%   Grids must follow definition from 'grid_structure'.  Both input grids
%   must be strictly two-dimensional both in terms of `griddim` and in
%   terms of size(nodes.coords, 2).
%
%   If the two grids do not share any common edges, the resulting combined
%   grid will represent a topologically disconnected grid.
%
% RETURNS:
%   G - Resultant combined grid from `G1` and `G2`.
%
% NOTE:
%   The result grid (`G`) does not provide derived geometric primitives
%   (e.g., cell volumes).  Such information must be explicitly computed
%   through a subsequent call to function `computeGeometry`.
%
% SEE ALSO:
%   `computeGeometry`.

%{
Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.

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

   opt = struct('tol', 1.0e-3);
   opt = merge_options(opt, varargin{:});

   % Insert extra nodes to ensure boundary conformity between the two grids
   G1_conformal = ensure_conformal_boundaries(G1, G2, opt.tol);
   G2_conformal = ensure_conformal_boundaries(G2, G1, opt.tol);

   % making combined grid (with redundant nodes and faces)
   G = concatenate_grids(G1_conformal, G2_conformal);

   % compute boundary loops (assuming no holes, for now)
   % a boundary loop is a list of oriented nodes and faces that constitutes
   % a grid's boundary.
   loop1 = boundary_loop(G1_conformal);
   loop2 = boundary_loop(G2_conformal);

   % identify coinciding nodes
   cnodes = coinciding_nodes(G1_conformal, loop1.nodes, ...
                             G2_conformal, loop2.nodes, opt.tol);

   % purging redundant nodes
   keepnodes = cnodes(:,1);
   discardnodes = cnodes(:,2) + G1_conformal.nodes.num;

   % faces
   issparse = true;
   m1 = accumarray([(1:numel(G.faces.nodes))', G.faces.nodes], 1, ...
                   [numel(G.faces.nodes), G.nodes.num], [], 0, issparse);
   subst = (1:G.nodes.num)';
   subst(discardnodes) = keepnodes;
   subst = accumarray([(1:G.nodes.num)', subst], 1, [], [], 0, issparse);
   subst(:, sum(subst) == 0) = [];

   [i, ~] = find((m1 * subst)');
   G.faces.nodes = i;

   % nodes
   G.nodes.num = G.nodes.num - numel(discardnodes);
   G.nodes.coords(discardnodes, :) = [];

   % Merge faces, i.e. ensuring neighbor cells share the same faces.
   % Result array 'mergefaces' is m-by-2.
   mergefaces = identify_mergefaces(loop1, loop2, cnodes, G1_conformal.faces.num);

   if ~isempty(mergefaces)
      G = merge_all_faces(G, mergefaces);
   end

   % Purging remaining, unused faces
   G = cleanup_unused_faces(G);
end

%--------------------------------------------------------------------------

function G = concatenate_grids(G1, G2)
   % make simple concatenated grid, without regards for duplication of
   % nodes/faces

   % nodes
   G.nodes.num = G1.nodes.num + G2.nodes.num;
   G.nodes.coords = [G1.nodes.coords; G2.nodes.coords];

   % faces
   G.faces.num = G1.faces.num + G2.faces.num;
   G.faces.nodePos = (1:G.faces.num + 1)' * 2 - 1;

   neigh1 = G1.faces.neighbors;
   neigh2 = G2.faces.neighbors;
   neigh2(neigh2 ~= 0) = neigh2(neigh2 ~=0) + G1.cells.num;

   G.faces.neighbors = [neigh1; neigh2];
   G.faces.nodes = [G1.faces.nodes; G2.faces.nodes + G1.nodes.num];

   % cells
   G.cells.num = G1.cells.num + G2.cells.num;
   G.cells.facePos = [G1.cells.facePos; G2.cells.facePos(2:end) + G1.cells.facePos(end)-1];
   G.cells.faces = [G1.cells.faces(:,1); G2.cells.faces(:,1) + G1.faces.num];
   G.cells.indexMap = []; % @@ unused, but added for backwards compatibility

   % other
   G.griddim = 2;
   G.type = { mfilename };
end

%--------------------------------------------------------------------------

function G = ensure_conformal_boundaries(targetgrid, othergrid, tol)
% loop through the boundary of 'targetgrid' and insert extra nodes (and
% split the resulting faces) wherever faces are intersected by boundary
% nodes from 'othergird' (away from its own boundary nodes).  The end
% result is a grid where the only boundary intersections with 'othergrid'
% are in terms of overlapping nodes, never nodes that intersect faces.

   boundary_faces = @(G) find(any(G.faces.neighbors == 0, 2));

   boundary_nodes0 = @(G, bf) ...
      unique(G.faces.nodes(mcolon(G.faces.nodePos(bf), ...
                                  G.faces.nodePos(bf + 1) - 1)));
   
   boundary_nodes = @(G) boundary_nodes0(G, boundary_faces(G));

   bfaces       = boundary_faces(targetgrid);
   bnodes       = boundary_nodes0(targetgrid, bfaces);
   bnodes_other = boundary_nodes(othergrid);

   % eliminate nodes that are already coinciding between the two grids
   cnodes = coinciding_nodes(othergrid, bnodes_other, targetgrid, bnodes, tol);
   bnodes_other = setdiff(bnodes_other, cnodes(:,1));

   bo_coords = othergrid.nodes.coords(bnodes_other, :);

   % intersections between nodes in 'bnodes_other' and faces in 'bfaces'
   [ifaces, ipoints] = point_face_isect(bfaces, targetgrid, bo_coords, tol);

   % inserting new nodes and faces
   G = targetgrid;
   G.cells.faces = G.cells.faces(:,1); % remove 'tags'
   unique_ifaces = unique(ifaces);

   for ui = reshape(unique_ifaces, 1, [])
      % determine point(s) to insert
      ins_pts = ipoints(ifaces == ui, :);
      num_new_pts = size(ins_pts, 1);

      % determine affected cell (should be only one), and sorting new
      % points along the affected edge
      assert(prod(G.faces.neighbors(ui, :), 2) == 0);
      cell = sum(G.faces.neighbors(ui, :));
      face_nodes = G.faces.nodes(mcolon(G.faces.nodePos(ui), G.faces.nodePos(ui+1)-1));
      ins_pts = order_points(ins_pts, G.nodes.coords(face_nodes(1), :));

      % updating nodes
      old_nodenum = G.nodes.num;
      G.nodes.num = G.nodes.num + num_new_pts;
      G.nodes.coords = [G.nodes.coords; ins_pts];

      % updating faces (add new faces - the old one will be purged later)
      old_facenum = G.faces.num;
      G.faces.num = G.faces.num + num_new_pts + 1;

      seq = repmat(old_nodenum+1:G.nodes.num, 2, 1);
      seq = [face_nodes(1); seq(:); face_nodes(2)];
      if G.faces.neighbors(ui, 1) == 0
         % face normal pointing into cell -> reverse order of replacement
         % faces to preserve consistent ordering of faces within cell
         seq = reshape(fliplr(reshape(seq, 2, [])),[] ,1);
      end

      G.faces.nodes = [G.faces.nodes; seq];
      G.faces.nodePos = [G.faces.nodePos; ...
                         (2:2:2*(num_new_pts+1))' + G.faces.nodePos(end)];
      G.faces.neighbors = [G.faces.neighbors; ...
                          repmat(G.faces.neighbors(ui,:), num_new_pts + 1, 1)];

      % updating cells
      facenums = diff(G.cells.facePos);
      f_ix = find(G.cells.faces == ui);  assert(numel(f_ix) == 1);
      G.cells.faces = [G.cells.faces(1:f_ix-1);
                       (old_facenum + 1 : old_facenum + num_new_pts + 1)';
                       G.cells.faces(f_ix+1:end)];
      facenums(cell) = facenums(cell) + num_new_pts;
      G.cells.facePos = cumsum([1; facenums(:)]);
   end

   G = cleanup_unused_faces(G);
end

%--------------------------------------------------------------------------

function pts = order_points(pts, ref_point)
% order points in increasing distance from reference point

   dists = sqrt(sum((pts - ref_point).^2, 2));
   [~, order] = sort(dists);
   pts = pts(order, :);
end

%--------------------------------------------------------------------------

function [ifaces, ipoints] = point_face_isect(faces, G, pts, tol)
% Detect intersections (within 'tol') between the indicated faces of G and
% the points whose 2D coords are given by 'pts'.

   nstart = G.faces.nodes(G.faces.nodePos(faces));
   nend   = G.faces.nodes(G.faces.nodePos(faces) + 1);

   s1 = G.nodes.coords(nstart,:);
   s2 = G.nodes.coords(nend, :);

   % unit vector along face direction
   fdir = s2 - s1;
   fdir = fdir ./ sqrt(sum(fdir.^2, 2)) ; % normalize

   [ifaces, ipoints] = deal([]);
   for p = pts'
      % vector from first face node to point (not normalized)
      pdir = p' - s1;

      % compute projected distance from 'p' to line passing through s1 and s2
      delta = fdir(:,2) .* pdir(:,1) - pdir(:,2) .* fdir(:, 1);

      ixs = find(abs(delta) < tol);
      if ~isempty(ixs)
         % point may be on a line on which several faces lie, but in the
         % end, only (zero or) one face should be intersected

         for i = reshape(ixs, 1, [])
            % check that intersection takes place within segment
            s = fdir(i,:);
            v1 = p' - s1(i,:);
            v2 = p' - s2(i,:);

            if (sum(v1 .* s) > 0) && (sum(v2 .* (-s)) > 0)
               ifaces = [ifaces; faces(i)]; %#ok
               ipoints = [ipoints; p']; %#ok
               break; % we found the intersection, there shouldn't be more
            end
         end
      end
   end
end

%--------------------------------------------------------------------------

function G = cleanup_unused_faces(G)
% Remove unused faces from G (i.e. those not referenced by any cell), and
% re-number the remaining faces to avoid gaps in numbering.

   missing_faces = setdiff(1:G.faces.num, G.cells.faces);

   % update G.faces
   G.faces.num = G.faces.num - numel(missing_faces);
   G.faces.neighbors(missing_faces, :) = [];

   num_nodes = diff(G.faces.nodePos);
   num_nodes(missing_faces) = [];
   G.faces.nodePos = cumsum([1; num_nodes]);
   nodes = reshape(G.faces.nodes, 2, [])';
   nodes(missing_faces, :) = [];
   nodes = nodes';
   G.faces.nodes = nodes(:);

   % Update G.cells
   N = max([missing_faces(:); G.cells.faces]);
   reindex = (1:N)';
   reindex(missing_faces) = [];
   new_ixs = nan(N,1);
   new_ixs(reindex) = 1:numel(reindex);
   G.cells.faces = new_ixs(G.cells.faces);
end

%--------------------------------------------------------------------------

function G = merge_all_faces(G, mfaces)
% Merge grid faces given in first column of 'mfaces' with the corresponding
% ones in the second column.  Leave both faces in place (although the ones
% in the first column will remain unused).  Unused faces can be purged
% later by a call to 'cleanup_unused_faces'.

   cells = sum(G.faces.neighbors(mfaces(:,1), :), 2);

   % Update faces.neighbors
   tmp = G.faces.neighbors(mfaces(:,2), :);
   tmp = reshape(tmp', [], 1);
   assert(numel(cells) == sum(tmp==0));
   tmp(tmp==0) = cells;
   tmp = reshape(tmp, 2, [])';
   G.faces.neighbors(mfaces(:,2), :) = tmp;

   % update G.cells.faces
   reindex = (1:max(G.cells.faces(:)))';
   reindex(mfaces(:, 1)) = mfaces(:,2);
   G.cells.faces = reindex(G.cells.faces);
end

%--------------------------------------------------------------------------

function mergefaces = identify_mergefaces(loop1, loop2, cnodes, offset)
% From the two loops 'loop1' and 'loop2', identify boundary face pairs that
% should be merged.  Such face-pairs are identified by each of the two
% faces in the pair having start- and end-nodes that coincide geometrically
% with that of the other (such nodes are already identified and provided by
% the 'cnodes' table).

   mergefaces = [];

   fnodes1 = [loop1.nodes, circshift(loop1.nodes, -1)];
   fnodes2 = [loop2.nodes, circshift(loop2.nodes, -1)];

   tally1 = fnodes1;
   tally2 = fnodes2;

   for c = cnodes'
      tally1(tally1 == c(1)) = 0;
      tally2(tally2 == c(2)) = 0;
   end

   f1_merges = find(sum(tally1, 2) == 0);
   f2_merges = find(sum(tally2, 2) == 0);
   assert(numel(f1_merges) == numel(f2_merges));

   for f1_cur = reshape(f1_merges, 1, [])
      % find counterpart face in loop2
      nodes1 = fnodes1(f1_cur,:);
      nodes2 = [cnodes(cnodes(:,1) == nodes1(2), 2), ...
                cnodes(cnodes(:,1) == nodes1(1), 2)];

      f2_cur = find(prod(fnodes2 == nodes2, 2));
      assert(numel(f2_cur) == 1); % there should be exactly one matching face

      mergefaces = [ mergefaces; ...
         [loop1.faces(f1_cur), loop2.faces(f2_cur) + offset] ...
      ]; %#ok<AGROW>
   end
end

%--------------------------------------------------------------------------

function cnodes = coinciding_nodes(G1, n1, G2, n2, tol)
% determine which of the nodes 'n1' in grid 'G1' and 'n2' in grid 'G2' that
% overlap within the given geometric tolerance.

   pts1 = G1.nodes.coords(n1, :);
   pts2 = G2.nodes.coords(n2, :);

   % @@ the below could likely be made more efficient
   dx2 = (pts1(:,1)' - pts2(:,1)).^2;
   dy2 = (pts1(:,2)' - pts2(:,2)).^2;

   dist = sqrt(dx2 + dy2);

   [i, j] = find(dist < tol);

   cnodes = [n1(j), n2(i)];
end

%--------------------------------------------------------------------------

function loop = boundary_loop(G)
% determine all boundary nodes and edges of a grid, and list them in
% counterclockwise order

   % identify all boundary faces
   bfaces = find(any(G.faces.neighbors == 0, 2));

   % all faces should have exactly two nodes
   assert(unique(diff(G.faces.nodePos)) == 2);

   bnodes = reshape(G.faces.nodes, 2, [])';
   bnodes = bnodes(bfaces,:);

   % loop as long as there are remaining faces
   n_indices = bfaces * 0;
   f_indices = bfaces * 0;

   f_indices(1) = bfaces(1);
   n_indices(1) = bnodes(1, 1);

   bfaces(1) = nan;

   count = 1;
   fix = 1;

   while true
      if n_indices(count) == bnodes(fix, 1)
         side = 2;
      else
         side = 1;
      end

      count = count + 1;
      n_indices(count) = bnodes(fix, side);

      bnodes(fix,:) = nan;

      fix = find(sum(bnodes == n_indices(count), 2), 1);

      f_indices(count) = bfaces(fix);
      bfaces(fix) = nan;
      if all(isnan(bfaces))
         break;
      end
   end

   if ~all(isnan(bfaces))
      % @@ support shouldn't be too hard to implement if needed
      error('Internal boundaries not yet supported');
   end

   % ensure clockwise orientation
   xy = G.nodes.coords(n_indices, :);

   if loop_orientation(xy(:,1), xy(:,2)) < 0
      n_indices = flipud(n_indices);
      f_indices = circshift(flipud(f_indices), -1);
   end

   loop.nodes = n_indices;
   loop.faces = f_indices;
end

%--------------------------------------------------------------------------

function orient = loop_orientation(x, y)
% Check if loop is clockwise (orient = 1) or counterclockwise (orient = -1)

   x = [x; x(1)];
   y = [y; y(1)];

   % ensure a pair number of elements
   if mod(numel(x), 2) == 1
      x = [x;x(1)];
      y = [y;y(1)];
   end

   i = (1:numel(x)-2)';
   i_p1 = i + 1;
   i_p2 = i + 2;

   area = sum (x(i_p1) .* (y(i_p2) - y(i)) + y(i_p1) .* (x(i) - x(i_p2))) / 2;

   orient = sign(area);
end
