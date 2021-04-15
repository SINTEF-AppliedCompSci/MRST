function varargout = plotPatches(G, faces, varargin)
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
%   faces   - Vector of face indices.  The graphical output of 'plotFaces'
%             will be restricted to the subset of grid faces from 'G'
%             represented by 'faces'.
%
%   colour  - Colour data specification.  Either a MATLAB 'ColorSpec'
%             (i.e., an RGB triplet (1-by-3 row vector) or a short or long
%             colour name such as 'r' or 'cyan'), or a PATCH
%             'FaceVertexCData' table suiteable for either indexed or
%             'true-colour' face colouring.  This data *MUST* be an m-by-1
%             column vector or an m-by-3 matrix.  We assume the following
%             conventions for the size of the colour data:
%
%                - ANY(SIZE(colour,1) == [1, NUMEL(faces)])
%                  One (constant) indexed colour for each face in 'faces'.
%                  This option supports 'flat' face shading only.  If
%                  SIZE(colour,1) == 1, then the same colour is used for
%                  all faces in 'faces'.
%
%                - SIZE(colour,1) == G.nodes.num
%                  One (constant) indexed colour for each node in 'faces'.
%                  This option must be chosen in order to support
%                  interpolated face shading.
%
%             OPTIONAL.  Default value: colour = 'y' (shading flat).
%
%   'pn'/pv - List of other property name/value pairs.  OPTIONAL.
%             This list will be passed directly on to function PATCH
%             meaning all properties supported by PATCH are valid.
%
% RETURNS:
%   h - Handle to resulting PATCH object.  The patch object is added to the
%       current AXES object.
%
% NOTES:
%   Function 'plotFaces' is implemented directly in terms of the low-level
%   function PATCH.  If a separate axes is needed for the graphical output,
%   callers should employ function newplot prior to calling 'plotFaces'.
%
% EXAMPLE:
%   % Plot grid with boundary faces on left side in red colour:
%   G     = cartGrid([5, 5, 2]);
%   faces = boundaryFaceIndices(G, 'LEFT', 1:5, 1:2, []);
%   plotGrid (G, 'faceColor', 'none'); view(3)
%   plotFaces(G, faces, 'r');
%
% SEE ALSO:
%   `plotCellData`, `plotGrid`, `newplot`, `patch`, `shading`.

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

   assert(~isCoarseGrid(G), ...
          'Do not send coarse grids to plotPatches');

   if isempty(faces)
      if nargout > 0, varargout(1 : nargout) = { [] }; end
      return
   end

   assert(isnumeric(faces) && all(isfinite(faces)) && all(faces > 0), ...
         ['Face list ''faces'' must be an array of ', ...
          'finite positive numbers.'])

   if mod(numel(varargin), 2) == 0
      colour   = 'yellow';
   else
      colour   = varargin{1};
      varargin = varargin(2 : end);
   end

   % Extract face topology for subset of grid faces and the coordinates for
   % the actual vertices ('verts') present in this topology.
   %
   dim = get_grid_dimension(G);

   if dim == 2
      nodes = getSortedCellNodes(G);
      pos   = G.cells.facePos;
   else
      nodes = G.faces.nodes;
      pos   = G.faces.nodePos;
   end

   [f, verts] = get_face_topo(nodes, pos, faces);
   v          = G.nodes.coords(verts, :);
   if isfield(G.nodes, 'z')
      v = [v, G.nodes.z(verts)];
   end

   if ~strcmpi(get(gcf, 'Renderer'), 'OpenGL') && ...
         need_render_fixes(varargin{:})
      set(gcf, 'Renderer', 'OpenGL');
   end

   % Massage colour data into form suitable for 'FaceVertexCData' property.
   %
   % From here, we assume that 'colour' is an m-by-1 (for indexed
   % colouring) or m-by-3 (for 'true colour' colouring) array.
   % Furthermore, 'm' is assumed to be either 1 (meaning all 'faces' should
   % be coloured using a single colour), NUMEL(faces) when the individual
   % faces are coloured separately, or G.nodes.num for interpolated face
   % colouring.
   %
   if ischar(colour), colour = get_rgb(colour); end

   if any(size(colour,1) == [1, numel(faces)])
      fc     = 'flat';
   elseif size(colour,1) == G.nodes.num
      colour = colour(verts,:);
      fc     = 'interp';
   end

   if size(f, 1) == size(v, 1)
      % plotPatches uses FaceVertexCData to plot face data.  MATLAB's PATCH
      % function decides if the data is given per vertex or face based on a
      % check of lengths.  However if the sizes are exactly *equal*, PATCH
      % defaults to interpreting the property as vertex data.  In MRST,
      % this can happen when dynamically plotting subsets of fronts.  Work
      % around this behaviour by adding a single dummy vertex with all NaN
      % coordinates.  This has no influence on the plot, but ensures that
      % the 'colour' data is correctly interpreted as face-related.
      v = [v; NaN([1, size(v, 2)])];
   end

   % Build final patch for graphical output (Note: added to GCA).
   if size(colour,1) == 1
      % Separate one-colour treatment to enable vector graphics (e.g.,
      % PRINT('-depsc2', ...)).
      h = patch('Faces'    , f     , 'Vertices', v, ...
                'FaceColor', colour, varargin{:});
   else
      h = patch('Faces'          , f      , 'Vertices' , v , ...
                'FaceVertexCData', colour , 'FaceColor', fc, varargin{:});
   end

   if nargout > 0, varargout{1} = h; end
end

%--------------------------------------------------------------------------
% Private helpers follow.
%--------------------------------------------------------------------------

function dim = get_grid_dimension(G)
   dim = 3;

   if all(diff(G.faces.nodePos) == 2) || ...
         (isfield(G, 'dim') && G.dim == 2)

      dim = 2;
   end
end

%--------------------------------------------------------------------------

function [f, present] = get_face_topo(nodes, pos, faces)
   eIX = pos;
   nn  = double(diff([pos(faces), ...
                      pos(faces + 1)], [], 2));
   fn  = double(nodes(mcolon(eIX(faces), eIX(faces + 1) - 1), 1));

   m   = numel(faces);
   n   = max(nn);
   f   = nan([n, m]);

   % Extract only those nodes/vertices actually present in the subset of
   % grid faces represented by 'faces'.  Create local numbering for these
   % vertices.
   %
   present           = false([max(nodes), 1]);
   present(fn)       = true;

   node_num          = zeros([max(nodes), 1]);
   node_num(present) = 1 : sum(double(present));

   off = reshape((0 : m - 1) .* n, [], 1);

   f(mcolon(off + 1, off + nn)) = node_num(fn);

   tmp         = isfinite(f);
   nnode       = sum(tmp,1);
   ind         = sub2ind(size(f),nnode,1:size(f,2));
   tmp         = repmat(f(ind),size(f,1),1);
   f(isnan(f)) = tmp(isnan(f));
   % PATCH requires that the 'Faces' property be a matrix of size
   % (number of faces)-by-(number of vertices).
   %
   f = f .';
end

%--------------------------------------------------------------------------

function b = need_render_fixes(varargin)
   b = false;

   if (ispc || ismac) && ~isempty(varargin)
      % Render fixes (possibly) needed on Microsoft Windows and MacOS X.
      opt = varargin(1 : 2 : end);
      val = varargin(2 : 2 : end);

      % Look for *LAST* occurrence of 'FaceAlpha' as that will determine
      % whether or not we need to override rendering options.
      %
      i = find(strncmpi('facea', opt, numel('facea')), 1, 'last');
      if ~isempty(i)
         % Render fixes needed for non-scalar 'FaceAlpha' or when
         % FaceAlpha < 1.
         b = (numel(val{i}) > 1) || ...
             ((numel(val{i}) == 1) && (val{i} < 1));
      end
   end
end

%--------------------------------------------------------------------------

function rgb = get_rgb(colour)
   switch lower(colour)
      case {'y', 'yellow' }, rgb = [1, 1, 0];
      case {'m', 'magenta'}, rgb = [1, 0, 1];
      case {'c', 'cyan'   }, rgb = [0, 1, 1];
      case {'r', 'red'    }, rgb = [1, 0, 0];
      case {'g', 'green'  }, rgb = [0, 1, 0];
      case {'b', 'blue'   }, rgb = [0, 0, 1];
      case {'w', 'white'  }, rgb = [1, 1, 1];
      case {'k', 'black'  }, rgb = [0, 0, 0];
      otherwise            , rgb = [0, 0, 1]; % Unknown colour -> 'blue'.
   end
end
