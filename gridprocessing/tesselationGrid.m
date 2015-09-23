function G = tesselationGrid(p, t)
%Construct valid grid definition from points and triangle list
%
% SYNOPSIS:
%   G = tesselationGrid(P, T)
%
% PARAMETERS:
%   P     - Node coordinates.  Must be an m-by-2 matrix, one row for each
%           node/point.
%
%   T     - Tesselation list: an n-by-k matrix where each row holds node
%           numbers for a k-polygon.
%
% RETURNS:
%   G     - Valid grid definition.
%
% EXAMPLE:
%
% SEE ALSO:
%   triangleGrid, tetrahedralGrid, grid_structure.

%{
Copyright 2009-2015 SINTEF ICT, Applied Mathematics.

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


   assert (size(p, 2) == 2, ...
          ['Function ''%s'' is only supported in two ', ...
           'space dimensions.'], mfilename);

   assert (nargin == 2,...
       'Must supply both a point set and a tesselation');

   assert (all(max(t) <= size(p, 1)), ...
       'Tesselation list ''T'' references invalid points.');

   assert (all(min(t) > 0), ...
       'Tesselation list ''T'' references invalid points.');

   n = size(t,2);
   idx               = rldecode([1:n 1]',[1 2*ones(1,n-1) 1]');
   [fn, i]           = sort(reshape(t(:, idx)', 2, []));
   [fn, cf, cf]      = unique(fn', 'rows');  %#ok

   G.faces.nodes     = reshape(fn', [], 1);
   G.cells.faces     = cf;

   G.nodes.coords    = p;
   G.nodes.num       = size(p, 1);

   G.cells.num       = size(t, 1);
   G.cells.facePos   = cumsum([1; repmat(n, [G.cells.num, 1])]);

   cellNo            = rldecode(1:G.cells.num, diff(G.cells.facePos), 2).';

   G.faces.num       = double(max(G.cells.faces));
   G.faces.neighbors = accumarray([G.cells.faces, i(1,:)'], cellNo, ...
                                  [G.faces.num, 2]);

   G.faces.nodePos   = cumsum([1; repmat(2, [G.faces.num, 1])]);

   G.type            = { mfilename };
   G.griddim         = 2;

   % Uniquify nodes
   h = zeros([G.nodes.num, 1]);
   h(G.faces.nodes) = 1;

   if sum(h) ~= G.nodes.num,
      [G.nodes.coords, a, map] = unique(G.nodes.coords, 'rows');  %#ok
      G.nodes.num = size(G.nodes.coords, 1);
      G.faces.nodes = map(G.faces.nodes);
   end
end
