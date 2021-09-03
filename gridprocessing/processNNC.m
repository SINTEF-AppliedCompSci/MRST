function nnc = processNNC(G, NNC, varargin)
%Establish explicit non-neighbouring connections from NNC keyword
%
% SYNOPSIS:
%   nnc = processNNC(G, NNC)
%   nnc = processNNC(G, NNC, excludeExisting)
%
% PARAMETERS:
%   G   - MRST grid as defined by `grid_structure`.
%
%   NNC - Explicit list of non-neighbouring (according to Cartesian cell
%         indices) as defined by function `readEclipseDeck` from ECLIPSE
%         keyword `NNC`.  Must be an m-by-7 array in which the first six
%         columns are (I,J,K) Cartesian indices of the first and second
%         connecting cells, respectively, and the seventh column is the
%         transmissibility of the corresponding connection.
%
%   excludeExisiting -
%         Flag indicating whether or not check if any of the explicit
%         non-neighbouring (tentatively non-geometrical) connections (pairs
%         of cells) are already present in the interface list and, if so,
%         exclude those from the new connections.
%
% RETURNS:
%   nnc - Structure describing explicit, additional non-neighbouring
%         connections.  Contains the following fields,
%
%           cells - Active cells connected across an NNC.  An m-by-2 array
%           of active cell numbers in the format of `G.faces.neighbors`.
%           Each row of nnc.cells represents a single non-neighbouring
%           connection.
%
%           trans - Transmissibility of the corresponding connections.
%
% NOTE:
%   Option `excludeExisting=true` is implemented in terms of function
%   `ismember(...,'rows')` which uses function `sortrows`.  This is
%   potentially very costly.
%
% SEE ALSO:
%   `readEclipseDeck`, `ismember`, `sortrows`.

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

   assert (isfield(G, 'cartDims') && (numel(G.cartDims) == 3), ...
           'Grid must be three-dimensional and logically Cartesian');

   assert (size(NNC, 2) == 7, ...
          ['Function ''%s'' only supports a simplified explicit ', ...
           'NNC list defintion.'], mfilename);

   s2i = @(ijk) sub2ind(G.cartDims, ijk(:,1), ijk(:,2), ijk(:,3));
   act = false([prod(G.cartDims), 1]);
   act(G.cells.indexMap) = true;

   i1 = s2i(NNC(:, 0*3 + (1:3)));
   i2 = s2i(NNC(:, 1*3 + (1:3)));
   p  = act(i1) & act(i2);

   nnc = struct('cells', [cart2active(G, i1(p)), ...
                          cart2active(G, i2(p))], ...
                'trans', NNC(p, 7));

   excludeExisting = false;
   if nargin > 2 && (numel(varargin{1}) == 1) && ...
         (islogical(varargin{1}) || isnumeric(varargin{1})),
      excludeExisting = varargin{1};
   end

   if excludeExisting,
      i = ~ismember(sort(nnc.cells, 2), ...
                    sort(double(G.faces.neighbors), 2), 'rows');

      nnc.cells = nnc.cells(i, :);
      nnc.trans = nnc.trans(i);

      if ~isempty(nnc.cells) && isfield(G, 'nnc'),
         o = [ G.nnc ];
         c = vertcat(o.cells);

         i = ~ismember(sort(nnc.cells, 2), sort(c, 2), 'rows');

         nnc.cells = nnc.cells(i, :);
         nnc.trans = nnc.trans(i);
      end
   end

   if isempty(nnc.cells),
      nnc = [];
   end
end
