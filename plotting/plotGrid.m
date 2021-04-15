function varargout = plotGrid(G, varargin)
% Plot exterior faces of grid to current axes.
%
% SYNOPSIS:
%       plotGrid(G)
%       plotGrid(G, 'pn1', pv1, ...)
%       plotGrid(G, cells)
%       plotGrid(G, cells, 'pn1', pv1, ...)
%   h = plotGrid(...)
%
% PARAMETERS:
%   G       - Grid data structure.
%
%   cells   - Vector of cell indices defining sub grid.  The graphical
%             output of function `plotGrid` will be restricted to the
%             subset of cells from `G` represented by `cells`.
%
%             If unspecified, function `plotGrid` will behave as if the
%             caller defined
%
%                 `cells = 1 : G.cells.num`
%
%             meaning graphical output will be produced for all cells in
%             the grid model `G`.  If `cells` is empty (i.e., if
%             `isempty(cells)`), then no graphical output will be produced.
%
% KEYWORD ARGUMENTS:
%
%   'Any'   - Additional keyword arguments will be passed directly on to
%             function `patch` meaning all properties supported by `patch`
%             are valid.
%
% RETURNS:
%   h  - Handle to resulting patch object.  The patch object is added
%        directly to the current `axes` object (`gca`).
%        OPTIONAL.  Only returned if specifically requested.  If
%        `isempty(cells)`, then `h==-1`.
%
% NOTES:
%   Function `plotGrid` is implemented directly in terms of the low-level
%   function `patch`.  If a separate axes is needed for the graphical output,
%   callers should employ function `newplot` prior to calling `plotGrid`.
%
% EXAMPLE:
%   G = cartGrid([10, 10, 5]);
%
%   % 1) Plot grid with yellow colour on faces (default):
%   figure, plotGrid(G, 'EdgeAlpha', 0.1); view(3)
%
%   % 2) Plot grid with no colour on faces (transparent faces):
%   figure, plotGrid(G, 'FaceColor', 'none'); view(3)
%
% SEE ALSO:
%   `plotCellData`, `plotFaces`, `patch`, `newplot`

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


if numel(G) > 1
   error(msgid('Grid:MultiComponent'), ...
         'Cannot plot more than one grid at a time.');
end

% Default to showing the entire grid model.
%
cells = (1 : G.cells.num) .';

if mod(numel(varargin), 2) == 1
   % Caller requested graphical output from a particular subset of the grid
   % cells.  Honour request, but treat an empty 'cells' argument as if the
   % caller requested the default behaviour (i.e., cells = 1:G.cells.num).
   %
   if isnumeric(varargin{1})
      cells = varargin{1};
   elseif islogical(varargin{1})
      cells = find(varargin{1});
   else
      error(['Second parameter ''cells'' must either be a list of ', ...
             'explicit cell indices or a logical mask into the '  , ...
             'grid''s cells.']);
   end

   % Strip 'cells' argument off of remaining input argument list.
   %
   varargin = varargin(2 : end);
end

if isempty(cells)
   warning(msgid('SubGrid:Empty'), ...
          'Empty cell selection in ''plotGrid''.  No graphics for you.')
   if nargout > 0, varargout{1} = -1; end
   return
end

% Assert that remaining arguments at least appear to be 'name'/value pairs
% intended for PATCH.
%
if ~isempty(varargin)
  assert (~mod(numel(varargin), 2) && ...
           all(cellfun(@ischar, varargin(1 : 2 : end))), ...
     'Additional arguments to plotGrid should be ''prop''/value pairs.');
end

if G.griddim == 3
   f = boundaryFaces(G, cells);
   subFaces = @getSubFaces;
else
   f = cells;
   subFaces = @getSubCells;
end

if G.griddim == 1
   plot(G.cells.centroids, ones(G.cells.num, 1), varargin{:});
elseif isCoarseGrid(G)
   % Separate otions: edge-related stuff is sent to plotFaceOutline
   ix           = rldecode(strncmpi('edge', varargin(1:2:end), 4), 2);
   outline_opts = varargin(ix);
   varargin     = varargin(~ix);
   ix           = rldecode(strncmpi('line', varargin(1:2:end), 4), 2);
   outline_opts = [outline_opts varargin(ix)];
   varargin     = varargin(~ix);

   h  = plotPatches(G.parent, subFaces(G, f), ...
      'EdgeColor', 'none', 'FaceColor', 'y', varargin{:});
   set(get(h, 'Parent'), 'ZDir', 'reverse');
   h  = [h; plotFaceOutline(G, f, outline_opts{:})];

else
   h = plotPatches(G, f, 'FaceColor', 'y', varargin{:});
   set(get(h, 'Parent'), 'ZDir', 'reverse');
end

if nargout > 0, varargout{1} = h; end
end

function subf = getSubFaces(G, f)
   ix   = mcolon(G.faces.connPos(f), G.faces.connPos(f+1)-1);
   subf = G.faces.fconn(ix);
end

function subcells = getSubCells(G, cells)
   ix        = false(G.cells.num, 1);
   ix(cells) = true;
   subcells  = find(ix(G.partition));
end
