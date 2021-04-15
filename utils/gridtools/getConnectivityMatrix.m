function A = getConnectivityMatrix(N, varargin)
%Derive global, undirected connectivity matrix from neighbourship relation.
%
% SYNOPSIS:
%   A = getConnectivityMatrix(N)
%   A = getConnectivityMatrix(N, incDiag)
%   A = getConnectivityMatrix(N, incDiag, nRows)
%
% PARAMETERS:
%   N       - Neighbourship relation as defined by function
%             `getNeighbourship`.  Must consist entirely of internal
%             connections (i.e., no boundary connections must be included).
%
%   incDiag - Flag to indicate whether or not to include the "self"
%             connections (i.e., the diagonal entries) in the connectivity
%             matrix.  LOGICAL.  Default value: `incDiag = false` (don't
%             include "self" connections, in which case the resulting
%             matrix is the same as the graph's adjacency matrix).
%
%   nRows   - Size (number of rows) of the resulting connectivity matrix.
%             Used to create the diagonal entries of the connectivity
%             matrix when `incDiag` is `true`.  Scalar integer.  Default
%             value: nRows = -1 (determine matrix size from maximum index
%             in the neighbourship relation, `N`).
%
% RETURNS:
%   A - Undirected connectivity matrix that represents all cell's
%       neighbourship relations (edges).
%
% EXAMPLE:
%   % Compute and plot the connectivity ("graph") matrix resulting from the
%   % (internal) geometrical connections in a Cartesian grid.  Include the
%   % diagonal entries of the connectivity matrix.
%   %
%   [incBdry, incDiag] = deal(false, true);
%   G = cartGrid([60, 220, 85]);
%   N = getNeighbourship(G, 'Geometrical', incBdry);
%   A = getConnectivityMatrix(N, incDiag, G.cells.num);
%   spy(A)
%
% SEE ALSO:
%   `getNeighbourship`.

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


   incDiag = false;
   if nargin > 1 && (numel(varargin{1}) == 1) && ...
         (isnumeric(varargin{1}) || islogical(varargin{1}))
      incDiag = varargin{1};
   end

   nRows = -1;
   if nargin > 2 && (numel(varargin{2}) == 1) && isnumeric(varargin{2})
      nRows = varargin{2};
   end

   % Implementation note: We use ACCUMARRAY with ISSPARSE=TRUE to
   % automatically cater to index types other than DOUBLE (such as when the
   % neighbourship is the 'G.faces.neighbors' field from a 'cartGrid').
   %
   subs = [ reshape(N, [], 1) , reshape(fliplr(N), [], 1) ];

   if incDiag
      if nRows > 0
         m = nRows;
      else
         m = max(subs(:,1));
      end

      i    = reshape(1 : m, [], 1);
      subs = [ subs ; [ i , i ] ];
   end

   A = accumarray(subs, 1, [], [], [], true);
end
