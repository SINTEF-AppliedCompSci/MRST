function varargout = plotFaces(G, faces, varargin)
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
%             G.faces.num.  The graphical output of 'plotFaces'
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
%             As a special case function plotFaces supports a separate
%             boolean option 'outline' which, when set, draws the outline
%             edge of the 'faces' input argument.  The outline is defined
%             as those edges that appear exactly once in the edge list
%             implied by 'faces'.
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
%   plotCellData, plotGrid, newplot, patch, shading.

%{
Copyright 2009-2014 SINTEF ICT, Applied Mathematics.

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


   if isempty(faces),
      if nargout > 0, varargout(1 : nargout) = { [] }; end
      return
   end

   if ~isnumeric(faces),
      if islogical(faces) && numel(faces) == G.faces.num
          faces = find(faces);
      else
          error(msgid('FaceList:NonNumeric'), ...
                'Face list ''faces'' is not numeric or a boolean per face.')
      end
   end

   if mod(numel(varargin), 2) == 0,
      colour   = 'yellow';
   else
      colour   = varargin{1};
      varargin = varargin(2 : end);
   end
   
   [plotOutline, varargin] = do_outline_p(varargin{:});


   % Extract face topology for subset of grid faces and the coordinates for
   % the actual vertices ('verts') present in this topology.
   %
   face_topo  = str2func(sprintf('get_face_topo_%dd', G.griddim));
   [f, verts] = face_topo     (G, faces);
   v          = G.nodes.coords(verts, :);

   if isfield(G.nodes, 'z'),
      assert (size(v,2) == 2, ...
             ['Vertex Z coordinate cannot be specified in ', ...
              'a general 3D grid']);

      v = [v, G.nodes.z(verts)];
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

   if any(size(colour,1) == [1, numel(faces)]),
      fc     = 'flat';
   elseif size(colour,1) == G.nodes.num,
      colour = colour(verts,:);
      fc     = 'interp';
   end

   if size(f, 1) == size(v, 1)
   % plotFaces uses FaceVertexCData to plot face data. The internal MATLAB
   % patch function decides if the data is given per vertex or face based
   % on a check of lengths. However, if the sizes are exactly *equal* it
   % defaults to vertex data. This can happen when dynamically plotting
   % subsets of fronts. In this case, we add a dummy vertex consisting of
   % NaN which has no influence on the plot, but ensures that the data is
   % interpreted as face data.
       v = [v; NaN(1, size(v, 2))];
   end

   % Build final patch for graphical output (Note: added to GCA).
   if size(colour,1) == 1 && (numel(faces) ~= 1 || size(colour,2) == 3)
      % Separate one-colour treatment to enable vector graphics (e.g.,
      % PRINT('-depsc2', ...)).

      h = patch('Faces'    , f     , 'Vertices', v, ...
                'FaceColor', colour, varargin{:});
   else
      h = patch('Faces'          , f      , 'Vertices' , v , ...
                'FaceVertexCData', colour , 'FaceColor', fc, varargin{:});
   end
   set(get(h, 'Parent'), 'ZDir', 'reverse')

   if plotOutline,
      pts = findFaceOutline(G, faces);
      do_hold = ishold();
      hold on, plot3(pts(:,1), pts(:,2), pts(:,3), 'k');
      if ~do_hold, hold off, end
   end

   if nargout > 0, varargout{1} = h; end
end

%--------------------------------------------------------------------------
% Private helpers follow.
%--------------------------------------------------------------------------

function [f, present] = get_face_topo_2d(G, cells)  %#ok
   eIX = G.cells.facePos;
   nn  = double(diff([G.cells.facePos(cells), ...
                      G.cells.facePos(cells + 1)], [], 2));

   cellNodes = getSortedCellNodes(G);

   cn  = double(cellNodes(mcolon(eIX(cells), eIX(cells + 1) - 1), 1));

   m   = numel(cells);
   n   = max(nn);
   f   = nan([n, m]);

   % Extract only those nodes/vertices actually present in the subset of
   % grid faces represented by 'faces'.  Create local numbering for these
   % vertices.
   %
   present           = false([G.nodes.num, 1]);
   present(cn)       = true;

   node_num          = zeros([G.nodes.num, 1]);
   node_num(present) = 1 : sum(double(present));

   off = reshape((0 : m - 1) .* n, [], 1);

   f(mcolon(off + 1, off + nn)) = node_num(cn);

   % PATCH requires that the 'Faces' property be a matrix of size
   % (number of faces)-by-(number of vertices).
   %
   f = f .';
end

%--------------------------------------------------------------------------

function [f, present] = get_face_topo_3d(G, faces)  %#ok
   eIX = G.faces.nodePos;
   nn  = double(diff([G.faces.nodePos(faces), ...
                      G.faces.nodePos(faces + 1)], [], 2));
   fn  = double(G.faces.nodes(mcolon(eIX(faces), eIX(faces + 1) - 1), 1));

   m   = numel(faces);
   n   = max(nn);
   f   = nan([n, m]);

   % Extract only those nodes/vertices actually present in the subset of
   % grid faces represented by 'faces'.  Create local numbering for these
   % vertices.
   %
   present           = false([G.nodes.num, 1]);
   present(fn)       = true;

   node_num          = zeros([G.nodes.num, 1]);
   node_num(present) = 1 : sum(double(present));

   off = reshape((0 : m - 1) .* n, [], 1);

   f(mcolon(off + 1, off + nn)) = node_num(fn);

   tmp=isfinite(f);
   nnode=sum(tmp,1);
   ind=sub2ind(size(f),nnode,1:size(f,2));
   tmp=repmat(f(ind),size(f,1),1);
   f(isnan(f))=tmp(isnan(f));
   % PATCH requires that the 'Faces' property be a matrix of size
   % (number of faces)-by-(number of vertices).
   %
   f = f .';
end

%--------------------------------------------------------------------------

function [plotOutline, varargin] = do_outline_p(varargin)
   % Does caller of 'plotFaces' request an outline plot?

   opt = struct('outline', false);
   [opt, varargin] = merge_options(opt, varargin{:});

   plotOutline = opt.outline;
end

%--------------------------------------------------------------------------

function rgb = get_rgb(colour)
   switch lower(colour),
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

   if numel(nodes) > 0,
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
