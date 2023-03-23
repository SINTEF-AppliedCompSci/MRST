function varargout = plotGrid_DFM(G, varargin)
%Plot exterior grid faces to current axes (reversed Z axis).
%
% This function has been modified from the original plotGrid to account for
% hybrid cells.
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
%             output of function 'plotGrid' will be restricted to the
%             subset of cells from 'G' represented by 'cells'.
%
%             If unspecified, function 'plotGrid' will behave as if the
%             caller defined
%
%                 cells = 1 : G.cells.num
%
%             meaning graphical output will be produced for all cells in
%             the grid model 'G'.  If 'cells' is empty (i.e., if
%             ISEMPTY(cells)), then no graphical output will be produced.
%
%   'pn'/pv - List of property names/property values.  OPTIONAL.
%             This list will be passed directly on to function PATCH
%             meaning all properties supported by PATCH are valid.
%
% RETURNS:
%   h  - Handle to resulting patch object.  The patch object is added
%        directly to the current AXES object (GCA).
%        OPTIONAL.  Only returned if specifically requested.  If
%        ISEMPTY(cells), then h==-1.
%
% NOTES:
%   Function 'plotGrid' is implemented directly in terms of the low-level
%   function PATCH.  If a separate axes is needed for the graphical output,
%   callers should employ function newplot prior to calling 'plotGrid'.
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
%   `plotCellData`, `plotFaces`, `patch`, `newplot`.

%{
Copyright 2009, 2010, 2011 SINTEF ICT, Applied Mathematics.

Portions Copyright 2011-2012 University of Bergen.

This file is part of DFM module of The MATLAB Reservoir Simulation Toolbox
(MRST).

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


if numel(G) > 1,
   error(msgid('Grid:MultiComponent'), ...
         'Cannot plot more than one grid at a time.');
end

% Default to showing the entire grid model.
%
cells = (1 : G.cells.num) .';

if mod(numel(varargin), 2) == 1,
   % Caller requested graphical output from a particular subset of the grid
   % cells.  Honour request, but treat an empty 'cells' argument as if the
   % caller requested the default behaviour (i.e., cells = 1:G.cells.num).
   %
   if isnumeric(varargin{1}),
      cells = varargin{1};
   end

   % Strip 'cells' argument off of remaining input argument list.
   %
   varargin = varargin(2 : end);
end

if isempty(cells),
   warning(msgid('SubGrid:Empty'), ...
          'Empty cell selection in ''plotGrid''.  No graphics for you.')
   if nargout > 0, varargout{1} = -1; end
   return
end

% Assert that remaining arguments at least appear to be 'name'/value pairs
% intended for PATCH.
%
if ~isempty(varargin)
  assert (all(cellfun(@ischar, varargin(1 : 2 : end))), ...
     'Additional arguments to plotGrid should be ''prop''/value pairs.');
end


if all(diff(G.faces.nodePos)==2),
   % For 2D grids, the faces to plot are the actual individual grid cells.
   f = cells;
elseif (size(G.nodes.coords,2)==2 && isfield(G.cells,'hybrid'))
   f = cells(G.cells.hybrid(cells)==0);
else
    f = boundaryFaces(G, cells);
end

h = plotFaces_DFM(G, f, 'EdgeColor', 'k', 'FaceColor', 'y', varargin{:});

if nargout > 0, varargout{1} = h; end
