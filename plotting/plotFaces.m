function varargout = plotFaces(G, varargin)
%Plot selection of coloured grid faces to current axes (reversed Z axis).
%
% SYNOPSIS:
%       plotFaces(G, faces)
%       plotFaces(G, faces, 'pn1', pv1, ...)
%       plotFaces(G, faces, colour)
%       plotFaces(G, faces, colour, 'pn1', pv1, ...)
%   h = plotFaces(...)
%
% PARAMETERS:
%   G       - Grid data structure.
%
%   faces   - Vector of face indices or a logical vector of length
%             G.faces.num.  The graphical output of `plotFaces`
%             will be restricted to the subset of grid faces from `G`
%             represented by `faces`.
%
%   colour  - Colour data specification.  Either a MATLAB `ColorSpec`
%             (i.e., an RGB triplet (1-by-3 row vector) or a short or long
%             colour name such as 'r' or 'cyan'), or a `patch`
%             `FaceVertexCData` table suitable for either indexed or
%             'true-colour' face colouring.  This data *MUST* be an m-by-1
%             column vector or an m-by-3 matrix.  We assume the following
%             conventions for the size of the colour data:
%
%                - `any(size(colour,1) == [1, numel(faces)])`
%                  One (constant) indexed colour for each face in `faces`.
%                  This option supports `flat` face shading only.  If
%                  `size(colour,1) == 1`, then the same colour is used for
%                  all faces in `faces`.
%
%                - `size(colour,1) == G.nodes.num`
%                  One (constant) indexed colour for each node in `faces`.
%                  This option must be chosen in order to support
%                  interpolated face shading.
%
%             OPTIONAL.  Default value: `colour = 'y'` (shading flat).
%
% KEYWORD ARGUMENTS:
%
%  'Any '    - Additional keyword arguments will be passed directly on to
%              function `patch` meaning all properties supported by `patch`
%              are valid.
%
%  'Outline' - Boolean option. When enabled, `plotFaces` draws the outline
%              edge of the `faces` input argument.  The outline is defined
%              as those edges that appear exactly once in the edge list
%              implied by `faces`.
%
% RETURNS:
%   h - Handle to resulting `patch` object.  The patch object is added to the
%       current `axes` object.
%
% NOTES:
%   Function `plotFaces` is implemented directly in terms of the low-level
%   function `patch`.  If a separate axes is needed for the graphical output,
%   callers should employ function `newplot` prior to calling `plotFaces`.
%
% EXAMPLE:
%   % Plot grid with boundary faces on left side in red colour:
%   G     = cartGrid([5, 5, 2]);
%   faces = boundaryFaceIndices(G, 'LEFT', 1:5, 1:2, []);
%   plotGrid (G, 'faceColor', 'none'); view(3)
%   plotFaces(G, faces, 'r');
%
% SEE ALSO:
%   `plotCellData`, `plotGrid`, `newplot`, `patch`, `shading`

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

   % Split input and options, assuming first string input is an option
   i     = cellfun(@ischar, varargin);
   first = find(i, 1, 'first');
   if mod(numel(varargin)-first+1, 2), first = first + 1; end

   if any(first)
      other_input = varargin(1:first-1);
      varargin    = varargin(first:end);
   else
      other_input = varargin;
      varargin = {};
   end

   % There should be an even number of elements in varargin here:
   assert(~mod(numel(varargin), 2), 'Huh?!');
   
   [plotOutline, varargin] = do_outline_p(varargin{:});
   
   color = 'y';
   switch numel(other_input)
      case 0
         if G.griddim == 2
            faces = 1:G.faces.num;
         else
            faces = boundaryFaces(G, 1:G.cells.num);
         end

         if ~any(strcmpi(varargin, 'FaceColor'))
            varargin = [varargin, {'FaceColor', 'y'}];
         end
      case 1
         faces = other_input{1};
         if islogical(faces)
            assert(numel(faces) == G.faces.num);
            faces = find(faces);
         end
         if ~any(strcmpi(varargin, 'FaceColor'))
            varargin = [varargin, {'FaceColor', 'y'}];
         end
      case 2
         faces = other_input{1};
         if islogical(faces)
            faces = find(faces);
         end
         color = other_input{2};
      otherwise
         error('What!?');
   end
   
   if isempty(faces)
      if nargout > 0, varargout(1 : nargout) = { [] }; end
      return
   end
   assert (min(faces) > 0, 'Cannot plot zero or negative face numbers');
   assert (max(faces) <= G.faces.num, ...
           'Faces in face list exceed number of faces in grid.');

   if G.griddim == 3

      h = plotFaces3D(G, faces, color, varargin{:});

   else

      h = plotFaces2D(G, faces, varargin{:});

   end
   
   if plotOutline
      pts = findFaceOutline(G, faces);
      do_hold = ishold();
      hold on, plot3(pts(:,1), pts(:,2), pts(:,3), 'k');
      if ~do_hold, hold off, end
   end
   
   if nargout > 0
       varargout{1} = h;
   end
end

%--------------------------------------------------------------------------

function h = plotFaces3D(G, faces, color, varargin)
   % If G is a coarsegrid, lookup finegrid data in parent
   if isCoarseGrid(G)
      h = plotFaces3DCoarseGrid(G, faces, varargin{:});
   else
      h = plotPatches(G, faces, color, 'EdgeColor', 'k', varargin{:});
      set(get(h, 'Parent'), 'ZDir', 'reverse')
   end
end

%--------------------------------------------------------------------------

function h = plotFaces2D(G, faces, varargin)
   % If G is a coarsegrid, lookup finegrid data in parent
   if isCoarseGrid(G)
      CG = G;
      f  = getSubFaces(G, faces);
      G  = G.parent;
   else
      f  = faces;
   end

   % Separate otions: marker-related stuff is sent to separate plotting
   ix          = rldecode(strncmpi('marker', varargin(1:2:end), 6), 2);
   marker_opts = varargin(ix);
   varargin    = varargin(~ix);

   % Plot fine grid edges specified by either coarse or fine grid.
   ix          = mcolon(G.faces.nodePos(f), G.faces.nodePos(f+1)-1);
   edges       = reshape(G.faces.nodes(ix), 2, []) .';
   h           = plotLineSegments(G, edges, varargin{:});
   set(get(h, 'Parent'), 'ZDir', 'reverse');

   % Find unique endpoints - if varargin contains 'marker*'
   if numel(marker_opts) > 0
      try
         faceno = rldecode(faces(:), (CG.faces.connPos(faces+1)-CG.faces.connPos(faces))*2, 1);
      catch %#ok
         faceno = rldecode(faces(:), 2);
      end

      e     = reshape(G.faces.nodes, 2, [])';
      e     = reshape(e(f,:)', [], 1);
      [E,n] = rlencode(sortrows([faceno, e]));
      nodes = E(n==1,2);

      holdstate = ishold;

      hold on;
      patch('vertices', G.nodes.coords(nodes,:), ...
            'faces'   , (1:numel(nodes)) .', marker_opts{:});

      if ~holdstate, hold off; end
   end
end

%--------------------------------------------------------------------------

function h = plotFaces3DCoarseGrid(G, faces, varargin)
   [f, fno] = getSubFaces(G, faces);

   if mod(numel(varargin), 2), varargin{1} = varargin{1}(fno); end

   % Marker-related otions are collected in marker_opts
   ix          = rldecode(strncmpi('marker', varargin(1:2:end), 6), 2);
   %marker_opts = varargin(ix);
   varargin    = varargin(~ix);

   % Edge-related options are collected in edge_opts
   ix = rldecode(strncmpi('edge', varargin(1:2:end), 4)', 2) | ...
        rldecode(strncmpi('line', varargin(1:2:end), 4)', 2);

   edge_opts = varargin(ix);
   varargin  = varargin(~ix);

   h = plotPatches(G.parent, f, 'EdgeColor', 'none', varargin{:});
   set(get(h, 'Parent'), 'ZDir', 'reverse')

   h = [h; plotFaceOutline(G, faces, edge_opts{:})];

%{
   if numel(marker_opts) > 0
      CG = G;
      G  = G.parent;

      ff = CG.faces.fconn;
      ffno = rldecode(1:CG.faces.num, diff(CG.faces.connPos), 2)';

      % Bug here:  A node is shared by two (or more) faces in 2D and three
      % (or more faces) in 3D.  This definition does not coincide with
      % "shared by two edges" hack above.

      % This check must be done per block: Check implementation of
      % cellNodes!

      [d, p] = copyRowsFromPackedData(G.faces.nodes, G.faces.nodePos, ff);
      fedges = [d, d(rot(p, 1))];
      faceno = rldecode(ffno, diff(p));

      tmp    = rlencode(sortrows([faceno, fedges(:,1); faceno, fedges(:,2)], [2,1]));

      [n,n]  = rlencode(tmp(:,2));
      N      = rldecode(n, n);
      nodes  = rlencode(tmp(N > 2, :));

      holdstate = ishold;
      hold on;

      %plot3(G.nodes.coords(nodes,1), G.nodes.coords(nodes, 2), G.nodes.coords(nodes, 3), 'linestyle','none', marker_opts{:});

      if ~holdstate, hold off; end

      %warning('Nodes in 3d is currently not supported');
   end
%}
end

%--------------------------------------------------------------------------

function h = plotLineSegments(G, e, varargin)
% Plot all line segments given by node pairs in each row in e.
   e = unique([e ; fliplr(e)], 'rows');

   if isfield(G.nodes, 'z')
      % Grid is 2d, but contains z coordinates in extra field.
      coord = [G.nodes.coords, G.nodes.z];
   else
      % Expected case.
      coord = G.nodes.coords;
   end

   h = patch('Vertices', coord, 'Faces', e, varargin{:});
end

%--------------------------------------------------------------------------

function [subf, fno] = getSubFaces(G, f)
   ix   = mcolon(G.faces.connPos(f), G.faces.connPos(f + 1) - 1);
   subf = G.faces.fconn(ix);
   fno  = rldecode(1 : numel(f), ...
                   G.faces.connPos(f + 1) - G.faces.connPos(f), 2) .';
end

%--------------------------------------------------------------------------

%{
function ix = rot(pos, offset)
   num    = diff(pos);
   offset = mod(offset, num); % net offset
   ix     = zeros(max(pos)-1, 1);

   ix(mcolon(pos(1:end-1), pos(1:end-1) + num - offset - 1)) = ...
      mcolon(pos(1:end-1) + offset, pos(2:end) - 1);

   ix(mcolon(pos(1:end-1) + num - offset, pos(2:end) - 1)) = ...
      mcolon(pos(1:end-1), pos(2:end) - 1 - num + offset);
end
%}

%--------------------------------------------------------------------------

%{
function [d, p] = copyRowsFromPackedData(d, p, rows)
   d = d(mcolon(p(rows), p(rows+1) - 1));
   p = cumsum([1 ; double(p(rows+1) - p(rows))]);
end
%}

%--------------------------------------------------------------------------

function [plotOutline, varargin] = do_outline_p(varargin)
   % Does caller of 'plotFaces' request an outline plot?

   opt = struct('outline', false);
   [opt, varargin] = merge_options(opt, varargin{:});

   plotOutline = opt.outline;
end

%--------------------------------------------------------------------------

function pts = findFaceOutline(g, faces)
   assert (size(g.nodes.coords, 2) == 3);
   % Find points on border of collection of grid faces.
   if numel(faces)==0, pts = zeros(0,3); return; end

   cellNo = rldecode(1:g.cells.num, diff(g.cells.facePos), 2) .';
   sgn    = 2*(cellNo == g.faces.neighbors(g.cells.faces(:,1), 1)) - 1;
   faceNo = rldecode(1:g.faces.num, diff(g.faces.nodePos), 2) .';

   fi     = false(g.faces.num, 1); fi(faces) = true;
   fn     = g.faces.nodes(fi(faceNo));

   nodeflag     = false([g.nodes.num, 1]);
   nodeflag(fn) = true;

   fe = faceEdges(g);
   fe(sgn(faceNo)<0, :) = fe(sgn(faceNo)<0, [2,1]);

   fe = fe (fi(faceNo),:);
   nodes = double(fe(all(nodeflag(fe), 2), :));

   if numel(nodes) > 0
      % Remove edges which appear more than once.  These edges are
      % considered internal in the collection of faces.  The remaining
      % edges are on the outer boundary.
      [nodes, n] = rlencode(sortrows(sort(nodes, 2)));
      nodes(n>1,:) = [];
   end

   pts = nan(size(nodes, 1)*3, 3);
   pts (1:3:end,:) = g.nodes.coords(nodes(:,1),:);
   pts (2:3:end,:) = g.nodes.coords(nodes(:,2),:);
end

%--------------------------------------------------------------------------

function fe = faceEdges(g)
   fe = [g.faces.nodes, g.faces.nodes([2:end,1])];
   fe(g.faces.nodePos(2:end)-1,2) = fe(g.faces.nodePos(1:end-1),1);
end
