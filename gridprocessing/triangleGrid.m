function G = triangleGrid(p, varargin)
%Construct valid grid definition from points and triangle list
%
% SYNOPSIS:
%   G = triangleGrid(P)
%   G = triangleGrid(P, T)
%
% PARAMETERS:
%   P     - Node coordinates.  Must be an m-by-2 matrix, one row for each
%           node/point.
%
%   T     - Triangle list: an n-by-3 matrix where each row holds node
%           numbers for a triangle.
%
%           OPTIONAL.  Default value:  T = DELAUNAY(P(:,1), P(:,2))
%
% RETURNS:
%   G     - Valid grid definition.
%
% EXAMPLE:
%   N     = 10;
%   N1    = 2*N;
%   N2    = 3*ceil(N/2)-2;
%   [X Y] = meshgrid(0:1:N1, 0:1:N2);
%   X     = sqrt(3) / 2 * X;
%   Y(:,2:2:end)=Y(:,2:2:end)+0.5;
%   p     = [X(:), Y(:)];
%   t     = delaunayn(p);
%   G     = triangleGrid(p, t);
%   cla, plotGrid(G)
%
% SEE ALSO:
%   `delaunay`, `tetrahedralGrid`, `grid_structure`

%{
Copyright 2009-2023 SINTEF Digital, Mathematics & Cybernetics.

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

   if nargin < 2,
      t = delaunay(p(:,1), p(:,2));

   elseif isnumeric(varargin{1}) && (size(varargin{1}, 2) == 3),
      t = varargin{1};

      assert (all(max(t) <= size(p, 1)), ...
              'Triangle list ''T'' references invalid points.');

      assert (all(min(t) > 0), ...
              'Triangle list ''T'' references invalid points.');

   else
      error(['Input parameter ''T'' does not appear to be a valid ', ...
             'triangulation of the points ''P''.']);
   end

   [fn, i]           = sort(reshape(t(:, [1,2, 2,3, 3,1])', 2, []));
   [fn, cf, cf]      = unique(fn', 'rows');  %#ok

   G.faces.nodes     = reshape(fn', [], 1);
   G.cells.faces     = cf;

   G.nodes.coords    = p;
   G.nodes.num       = size(p, 1);

   G.cells.num       = size(t, 1);
   G.cells.facePos   = cumsum([1; repmat(3, [G.cells.num, 1])]);

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
