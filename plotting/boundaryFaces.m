function [f, varargout] = boundaryFaces(g, varargin)
%Extract boundary faces from set of grid cells.
%
% SYNOPSIS:
%    f     = boundaryFaces(G)
%    f     = boundaryFaces(G, cells)
%   [f, c] = boundaryFaces(...)
%
% PARAMETERS:
%   G     - Grid data structure.
%
%   cells - Non-empty subset of cells from which to extract boundary faces.
%           OPTIONAL.  Default value: `cells = 1 : G.cells.num`, meaning
%           all external faces for all grid cells will be extracted. This
%           amounts to extracting the entire boundary of 'G'.
%
% RETURNS:
%   f - List of faces bounding the sub domain given by `cells`.
%
%   c - List of specific grid cells connected to the individual faces in
%       `f`.  This may be useful for plotting cell data (e.g., the cell
%       pressure) on the sub domain faces by means of function `plotFaces`.
%
% EXAMPLE:
%   G    = cartGrid([40, 40, 5]);
%   rock = <some rock data structure for G>;
%
%   % 1) Plot (external) geometry of 'G'.
%   f  = boundaryFaces(G);
%   hg = plotFaces(G, f);
%
%   % 2) Plot horizontal permeability along diagonal of reservoir
%   [f, c] = boundaryFaces(G, 1 : G.cartDims(1) + 1 : G.cells.num);
%   hd     = plotFaces(G, f, log10(rock.perm(c,1)), 'FaceAlpha', 0.625);
%
% SEE ALSO:
%   `plotFaces`

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


   if (nargin > 1) && isnumeric(varargin{1})
      sub = varargin{1};
      if any(sub == 0)
         warning(msgid('Outside:InCellSubset'), ...
            'Cell zero (outside) included in subset. Ignored.');
      end
   elseif (nargin>1) && islogical(varargin{1})
      sub = find(varargin{1});
   else
       sub = (1 : g.cells.num) .';
   end

   present          = false(1, g.cells.num + 1);
   present(sub + 1) = true;
   present( 0  + 1) = false;  % Explicitly ignore cell '0' (outside).

   % Extract boundary faces:
   %   Those faces for which one of the connecting cells is 'present'
   %   (i.e., within the sub domain 'sub') and the other is ~present
   %   (outside the sub domain).  This means that
   %
   %       SUM(present(g.faces.neighbors), 2) == 1
   %
   %   for the external faces.  Add one to account for the outside being
   %   'cell 0'.
   %
   tmp = present(g.faces.neighbors + 1);
   f = find(xor(tmp(:,1), tmp(:,2)));

   if nargout > 1
      % User requested list of cells connected to the faces in 'f'.

      % We construct a column vector 'c' such that c(i) is the cell within
      % 'sub' which connects to the face f(i) for all i=1:numel(f).
      %
      % Algorithm:
      %   1) Construct reduced connection matrix 'n', one row for each face
      %      in 'f', such that one of the cells n(i,1) or n(i,2) is the
      %      required cell c(i).
      %
      %   2) Mask reduced connection matrix by cell presence.  Result is a
      %      matrix whose entries are either zero (where the original
      %      entries of 'n' contained a cell outside 'sub') or the original
      %      cell number (if inside 'sub').
      %
      %   3) Sum the rows of this matrix to obtain the corresponding cell
      %      number regardless of original column number.

      n = double(g.faces.neighbors(f,:));
      c = sum(n .* double(present(n + 1)), 2);

      assert (all(c > 0));

      varargout{1} = c;
   end
end
