function varargout = plotCellData(G, data, varargin)
%Plot exterior grid faces, coloured by given data, to current axes.
%
% SYNOPSIS:
%       plotCellData(G, data)
%       plotCellData(G, data, 'pn1', pv1, ...)
%       plotCellData(G, data, cells)
%       plotCellData(G, data, cells, 'pn1', pv1, ...)
%   h = plotCellData(...)
%
% PARAMETERS:
%   G       - Grid data structure.
%
%   data    - Scalar cell data with which to colour the grid.  One scalar,
%             indexed colour value for each cell in the grid or one
%             TrueColor value for each cell.  If a cell subset is specified
%             in terms of the 'cells' parameter, 'data' must either contain
%             one scalar value for each cell in the model or one scalar
%             value for each cell in this subset.
%
%   cells   - Vector of cell indices defining sub grid.  The graphical
%             output of function 'plotCellData' will be restricted to the
%             subset of cells from 'G' represented by 'cells'.
%
%             If unspecified, function 'plotCellData' will behave as if the
%             user defined
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
%   Function 'plotCellData' is implemented directly in terms of the
%   low-level function PATCH.  If a separate axes is needed for the
%   graphical output, callers should employ function newplot prior to
%   calling 'plotCellData'.
%
% EXAMPLE:
%   Given a grid 'G' and a reservoir solution structure 'resSol' returned
%   from, e.g., function 'solveIncompFlow', plot the cell pressure in bar:
%
%      figure, plotCellData(G, convertTo(resSol.pressure, barsa()));
%
% SEE ALSO:
%   plotFaces, boundaryFaces, patch, newplot.

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


if numel(G) > 1,
   error(msgid('Grid:MultiComponent'), ...
         'Cannot plot more than one grid at a time.');
end

assert (sum(size(data, 2) == [1, 3]) == 1, ...
       'Second input, DATA, must have one or three columns.');

% Default to providing graphical output from all cells in the grid model.
%
cells = (1 : G.cells.num) .';

if mod(numel(varargin), 2) == 1,
   % Caller requested graphical output from a particular subset of the grid
   % cells.  Honour that request, but only if it makes sense in context.
   %
   if isnumeric(varargin{1}),
      cells = varargin{1};
   elseif islogical(varargin{1}) && ...
         numel(varargin{1}) == G.cells.num,
      cells = find(varargin{1});
   else
      error(['Third parameter ''cells'' must either be a list of ', ...
             'explicit cell indices or a logical mask into the '  , ...
             'grid''s cells.']);
   end

   % Strip 'cells' argument off of remaining input argument list.
   %
   varargin = varargin(2 : end);
end

if isempty(cells),
   warning(msgid('SubGrid:Empty'), ...
          ['Empty cell selection in ''plotCellData''.', ...
           '  No graphics for you.'])
   if nargout > 0, varargout{1} = -1; end
   return
end

% Assert that remaining arguments at least appear to be 'name'/value pairs
% intended for PATCH.
%
if ~isempty(varargin),
   assert (iscellstr(varargin(1 : 2 : end)), ...
           'Additional parameters should be ''key''/value pairs.');
end

assert (ndims(data) == 2, ...
       ['Function ''%s'' only supports 2D data arrays ', ...
        '(i.e., vectors or matrices).'], mfilename);

assert (any(size(data, 2) == [1, 3]), ...
        'Function ''%s'' colour data must be m-by-1 or m-by-3 array.', ...
        mfilename);

assert (any(size(data, 1) == [G.cells.num, numel(cells)]),  ...
        'The DATA should have one value for each grid cell in output.');

if G.griddim == 3,
   % Define the boundary to be the boundary faces of the selected subset
   % where the values are not nan.
   selectcells = cells;
   datanan = any(isnan(data), 2);
   if any(datanan)
       selectcells = false(G.cells.num, 1);
       selectcells(cells) = true;
       selectcells = find(selectcells & ~datanan);
   end
   [f, c] = boundaryFaces(G, selectcells);
else
   % For 2D grids, the faces to plot are the actual individual grid cells.
   [f, c] = deal(cells);
end

if size(data, 1) < G.cells.num,
   renum        = zeros([G.cells.num, 1]);
   renum(cells) = 1 : numel(cells);
   c            = renum(c);

   assert (all(c > 0) && all(c <= size(data, 1)));
end

h = plotFaces(G, f, data(c, :), 'EdgeColor', 'none', varargin{:});

if nargout > 0, varargout{1} = h; end
