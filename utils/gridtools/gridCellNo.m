function cellno = gridCellNo(G, varargin)
%Construct map from half-faces to cells or cell subset
%
% SYNOPSIS:
%   cellno = gridCellNo(G)
%   cellno = gridCellNo(G, c)
%
% PARAMETERS:
%   G - Grid structure.
%
%   c - Cells for which to construct mapping.  Array of numeric cell
%       indices in the range 1:G.cells.num .
%
%       OPTIONAL.  If unspecified, function `gridCellNo` will behave as
%       if ::
%
%           c = 1 : G.cells.num
%
%       In other words the map will be constructed for all grid cells.
%
% RETURNS:
%   cellno -
%       Map from half-faces to cell or cell subset indices.  The
%       expression::
%
%           c(cellno(hf))
%
%       derives the cell in `G` to which half-face `hf` belongs.  This, in
%       turn, means that if `c` is unspecified, then `cellno(hf)` is that
%       cell directly.  The half-face `hf` must be relative to the cell
%       subset identified by `c` or global if `c` is not specified.

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

   mrstNargInCheck(1, 2, nargin);

   if nargin == 1,

      % Simple case.  Create 'cellno' for all cells.
      cellno = rldecode(1 : G.cells.num, diff(G.cells.facePos), 2) .';

   else

      % Caller specified cell subset.
      cells = varargin{1};

      assert (isnumeric(cells) && ~isempty(cells) && ...
              (1 <= min(cells)) && (max(cells) <= G.cells.num), ...
             ['Cell subset ''c'' must be numeric indices in the ', ...
              'range [1 .. %d]'], G.cells.num);

      nf = G.cells.facePos(cells + 1) - G.cells.facePos(cells);

      % Find the cell index of each face.
      cellno = rldecode(1 : numel(cells), nf, 2) .';

   end
end
